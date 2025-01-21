#include "logging.h"
#include "macros.h"

#include <stdbool.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/utsname.h>
#include <sys/resource.h>
#ifdef USE_MPI
#include <mpi.h>
#endif

#ifndef LOG_COMMIT
#define LOG_COMMIT "unknown"
#endif

#define LOG_STDOUT stdout
#define LOG_STDERR stdout

/* if you would like to use stderr output stream for errors, redefine above macro */

#ifdef LOG_COLOR
#define LOG_ALL_RESET     "\x1b[0m"
#define LOG_FGCOLOR_RESET "\x1b[39m"
#define LOG_BGCOLOR_RESET "\x1b[49m"
#define LOG_COLOR_ULINE   "\x1b[4m"
#define LOG_ULINE_RESET   "\x1b[24m"
#define LOG_COLOR_WHITE   "\x1b[1;37m"
#define LOG_COLOR_REDBG   "\x1b[1;97;41m"

static const char *level_colors[] =
    { LOG_COLOR_REDBG, LOG_COLOR_RED, LOG_COLOR_PURPLE, "", LOG_COLOR_GREEN, LOG_COLOR_BLUE };
#endif                          // LOG_COLOR

static const char *level_strings[] = { "FATAL", "ERROR", "WARN", "", "INFO", "DEBUG" };

static log_Level log_minlevel = LOG_INFO;
static FILE *fpstd = NULL;
static FILE *fperr = NULL;
static struct timespec ts_start;

int log_mpid_ = 0;

const char *log_get_date_()
{
    static char datbuf[128];
    time_t now = time(NULL);
    strftime(datbuf, 127, "%a %d %b %Y, %T %Z", localtime(&now));
    return (datbuf);
}

const char *log_get_nodename()
{
    static struct utsname kern;
    uname(&kern);
    return kern.nodename;
}

void log_init(FILE * fp)
{
#ifdef USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &log_mpid_);
#endif
    if (fp) {
        fpstd = fp;
        fperr = fp;
    } else {
        fpstd = LOG_STDOUT;
        fperr = LOG_STDERR;
    }
    clock_gettime(CLOCK_REALTIME, &ts_start);
}

void log_set_output(FILE * fp)
{
    if (fp) {
        fpstd = fp;
        fperr = fp;
    } else {
        fpstd = LOG_STDOUT;
        fperr = LOG_STDERR;
    }
}

void log_set_level(log_Level level)
{
    if (level < LOG_WARN)
        log_minlevel = LOG_ERROR;
    else
        log_minlevel = level;
}

void log_set_level_from_string(const char *s_level)
{
    bool b_match = false;

    /* SILENT: only output error, fatal */
    if (STRSTRCOMP(s_level, "SIL")) {
        log_warnc(0, "Switching to silent mode for logging. Only errors are reported from here on.\n");
        b_match = true;
        log_set_level(LOG_ERROR);
    }

    /* WARNING: only output warning, error, fatal */
    if (STRSTRCOMP(s_level, "WARN")) {
        b_match = true;
        log_set_level(LOG_WARN);
    }

    /* INFO: output info, warning, error, fatal */
    if (STRSTRCOMP(s_level, "INFO")) {
        b_match = true;
        log_set_level(LOG_INFO);
    }

    /* DEBUG: output everything */
    if (STRSTRCOMP(s_level, "DEB")) {
        b_match = true;
        log_set_level(LOG_DEBUG);
    }

    if (!b_match) {
        log_set_level(LOG_INFO);
        log_warn("Could not determine logging level from string '%s' - set to 'INFO' (default).", s_level);
    }
}

log_Level log_get_level()
{
    return log_minlevel;
}

const char *log_get_level_string(log_Level level)
{
    return level_strings[level];
}

int log_get_mpid()
{
    return log_mpid_;
}

void log_general_(log_Level level, const char *file, int line, const char *fmt, ...)
{
    if (level > log_minlevel)
        return;

    va_list ap;
    va_start(ap, fmt);

    time_t t;
    char buf[16];

    switch (level) {
      case LOG_DEBUG:
          t = time(NULL);
          buf[strftime(buf, sizeof(buf), "%H:%M:%S", localtime(&t))] = '\0';
#ifdef LOG_COLOR
          fprintf(fpstd, "%s%-5s%s[%s,%d] %s%s(%d)%s: ", level_colors[level], level_strings[level],
                  LOG_ALL_RESET, buf, log_mpid_, LOG_COLOR_BOLD, file, line, LOG_ALL_RESET);
#else
          fprintf(fpstd, "%-5s[%s,%d] %s(%d): ", level_strings[level], buf, log_mpid_, file, line);
#endif
          vfprintf(fpstd, fmt, ap);
          fflush(fpstd);
          break;
      case LOG_STD:
          vfprintf(fpstd, fmt, ap);
          fflush(fpstd);
          break;
      case LOG_INFO:
      case LOG_WARN:
          t = time(NULL);
          buf[strftime(buf, sizeof(buf), "%H:%M:%S", localtime(&t))] = '\0';
#ifdef LOG_COLOR
          fprintf(fpstd, "%s%-5s%s[%s,%d] ", level_colors[level], level_strings[level], LOG_ALL_RESET, buf, log_mpid_);
#else
          fprintf(fpstd, "%-5s[%s,%d] ", level_strings[level], buf, log_mpid_);
#endif
          vfprintf(fpstd, fmt, ap);
          fflush(fpstd);
          break;
      case LOG_ERROR:
          t = time(NULL);
          buf[strftime(buf, sizeof(buf), "%H:%M:%S", localtime(&t))] = '\0';
#ifdef LOG_COLOR
          fprintf(fperr, "%s%-5s%s[%s,%d] ", level_colors[level], level_strings[level], LOG_ALL_RESET, buf, log_mpid_);
#else
          fprintf(fperr, "%-5s[%s,%d] ", level_strings[level], buf, log_mpid_);
#endif
          vfprintf(fperr, fmt, ap);
          fflush(fperr);
          break;
      case LOG_FATAL:
          t = time(NULL);
          buf[strftime(buf, sizeof(buf), "%H:%M:%S", localtime(&t))] = '\0';
#ifdef LOG_COLOR
          fprintf(fperr, "%s%-5s[%s,%d]%s ", level_colors[level], level_strings[level], buf, log_mpid_, LOG_ALL_RESET);
#else
          fprintf(fperr, "%-5s[%s,%d] ", level_strings[level], buf, log_mpid_);
#endif
          vfprintf(fperr, fmt, ap);
          fflush(fperr);
#ifdef USE_MPI
          MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#else
          exit(EXIT_FAILURE);
#endif
          break;
    }

    va_end(ap);
}

void log_banner(log_Program prog)
{
    static struct utsname kern;
    uname(&kern);

    if (0 == log_mpid_) {
        log_std("\n***********************************************************************\n");
#ifdef LOG_COLOR
        if (prog == LOG_SOFI) {
            log_std("* %sSOFI: Parallel Viscoelastic Anisotropic Finite-Difference Modelling%s *\n",
                    LOG_COLOR_BOLD, LOG_ALL_RESET);
        } else if (prog == LOG_SNAP) {
            log_std("* %sSNAPMERGE: Merge of SOFI wavefield snapshot files                  %s *\n",
                    LOG_COLOR_BOLD, LOG_ALL_RESET);
        }
#else
        if (prog == LOG_SOFI) {
            log_std("* SOFI: Parallel Viscoelastic Anisotropic Finite-Difference Modelling *\n");
        } else if (prog == LOG_SNAP) {
            log_std("* SNAPMERGE: Merge of SOFI wavefield snapshot files                   *\n");
        }
#endif
        log_std("* written by Thomas Bohlen and colleagues                             *\n");
        log_std("* Geophysical Institute, KIT-Department of Physics                    *\n");
        log_std("* Karlsruhe Institute of Technology (KIT), Karlsruhe, Germany         *\n");
        log_std("* https://gpi.kit.edu/                                                *\n");
        log_std("* Program compiled from git commit: %-33s *\n", LOG_COMMIT);
        log_std("***********************************************************************\n");
    }
#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    log_info("Node: %s (%s %s, %s)\n", kern.nodename, kern.sysname, kern.machine, kern.release);

#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

#ifdef LOG_COLOR
    log_infoc(0, "%sProcess started %s.%s\n", LOG_COLOR_BOLD, log_get_date_(), LOG_ALL_RESET);
#else
    log_infoc(0, "Process started %s.\n", log_get_date_());
#endif
}

void log_finalize()
{
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);

    struct timespec ts_end;
    clock_gettime(CLOCK_REALTIME, &ts_end);

#ifdef LOG_COLOR
    log_infoc(0, "%sProcess terminated %s.%s\n", LOG_COLOR_BOLD, log_get_date_(), LOG_ALL_RESET);
#else
    log_infoc(0, "Process terminated %s.\n", log_get_date_());
#endif

#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    double rtime = (double)(ts_end.tv_sec - ts_start.tv_sec) + (double)(ts_end.tv_nsec - ts_start.tv_nsec) * 1.0e-9;
    double rhours = floor(rtime / 3600.0);
    double rrem_secs = rtime - rhours * 3600.0;
    double rmins = floor(rrem_secs / 60.0);
    double rsecs = rrem_secs - rmins * 60.0;
    char s_rhours[64], s_rmins[64], s_rsecs[64];
    if (rhours > 0.0)
        sprintf(s_rhours, "%ldh ", (long)rhours);
    else
        s_rhours[0] = '\0';
    if (rmins > 0.0 || rhours > 0.0)
        sprintf(s_rmins, "%ldm ", (long)rmins);
    else
        s_rmins[0] = '\0';
    sprintf(s_rsecs, "%lfs, user time: ", rsecs);
    double utime = (double)usage.ru_utime.tv_sec + (double)usage.ru_utime.tv_usec * 1.0e-6;
    double uhours = floor(utime / 3600.0);
    double urem_secs = utime - uhours * 3600.0;
    double umins = floor(urem_secs / 60.0);
    double usecs = urem_secs - umins * 60.0;
    char s_uhours[64], s_umins[64], s_usecs[64];
    if (uhours > 0.0)
        sprintf(s_uhours, "%ldh ", (long)uhours);
    else
        s_uhours[0] = '\0';
    if (umins > 0.0 || uhours > 0.0)
        sprintf(s_umins, "%ldm ", (long)umins);
    else
        s_umins[0] = '\0';
    sprintf(s_usecs, "%lfs, system time: ", usecs);
    double stime = (double)usage.ru_stime.tv_sec + (double)usage.ru_stime.tv_usec * 1.0e-6;
    double shours = floor(stime / 3600.0);
    double srem_secs = stime - shours * 3600.0;
    double smins = floor(srem_secs / 60.0);
    double ssecs = srem_secs - smins * 60.0;
    char s_shours[64], s_smins[64], s_ssecs[64];
    if (shours > 0.0)
        sprintf(s_shours, "%ldh ", (long)shours);
    else
        s_shours[0] = '\0';
    if (smins > 0.0 || shours > 0.0)
        sprintf(s_smins, "%ldm ", (long)smins);
    else
        s_smins[0] = '\0';
    sprintf(s_ssecs, "%lfs", ssecs);

    log_info("Real time: %s%s%s%s%s%s%s%s%s\n",
             s_rhours, s_rmins, s_rsecs, s_uhours, s_umins, s_usecs, s_shours, s_smins, s_ssecs);
}
