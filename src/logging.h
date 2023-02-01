
/** \file
 *  Functions to provide formatted logging to stdout/stderr or files.
 */

#ifndef __LOGGING_H__
#define __LOGGING_H__

#include <stdbool.h>
#include <stdio.h>

#ifdef LOG_COLOR
#define LOG_ALL_RESET     "\x1b[0m"
#define LOG_COLOR_BOLD    "\x1b[1m"
#define LOG_BOLD_RESET    "\x1b[21m"
#define LOG_COLOR_RED     "\x1b[1;31m"
#define LOG_COLOR_GREEN   "\x1b[1;32m"
#define LOG_COLOR_YELLOW  "\x1b[1;33m"
#define LOG_COLOR_BLUE    "\x1b[1;34m"
#define LOG_COLOR_PURPLE  "\x1b[1;35m"
#define LOG_COLOR_CYAN    "\x1b[1;36m"
#else
#define LOG_ALL_RESET     ""
#define LOG_COLOR_BOLD    ""
#define LOG_BOLD_RESET    ""
#define LOG_COLOR_RED     ""
#define LOG_COLOR_GREEN   ""
#define LOG_COLOR_YELLOW  ""
#define LOG_COLOR_BLUE    ""
#define LOG_COLOR_PURPLE  ""
#define LOG_COLOR_CYAN    ""
#endif

/** 
 * Enum to describe logging levels.
 */
typedef enum log_Level {
    LOG_FATAL = 0,
    LOG_ERROR,
    LOG_WARN,
    LOG_STD,
    LOG_INFO,
    LOG_DEBUG
} log_Level;

/** 
 * Enum to describe programs.
 */
typedef enum log_Program {
    LOG_SOFI = 0,
    LOG_SNAP = 1
} log_Program;

/* Output debug information (prefix DEBUG) */
#define log_debug(...) log_general_(LOG_DEBUG, __FILE__, __LINE__, ##__VA_ARGS__)

/* Output standard information (prefix INFO) */
#define log_info(...)  log_general_(LOG_INFO,  __FILE__, __LINE__, ##__VA_ARGS__)

/* Output standard information without prefix */
#define log_std(...)   log_general_(LOG_STD,   __FILE__, __LINE__, ##__VA_ARGS__)

/* Output warnings (prefix WARN) */
#define log_warn(...)  log_general_(LOG_WARN,  __FILE__, __LINE__, ##__VA_ARGS__)

/* Output errors (prefix ERROR) */
#define log_error(...) log_general_(LOG_ERROR, __FILE__, __LINE__, ##__VA_ARGS__)

/* Output fatal errors (prefix FATAL) and abort/shut down MPI processes */
#define log_fatal(...) log_general_(LOG_FATAL, __FILE__, __LINE__, ##__VA_ARGS__)

/* declaration only required to make the following macros 
 * work; this variable should never be used directly! */
extern int log_mpid_;

/* Output debug information (prefix DEBUG) only on given RANK */
#define log_debugc(RANK,...) \
  do { if((RANK)==log_mpid_) log_general_(LOG_DEBUG,  __FILE__, __LINE__, ##__VA_ARGS__); } while (0)

/* Output standard information (prefix INFO) only on given RANK */
#define log_infoc(RANK,...) \
  do { if((RANK)==log_mpid_) log_general_(LOG_INFO,  __FILE__, __LINE__, ##__VA_ARGS__); } while (0)

/* Output standard information without prefix only on given RANK */
#define log_stdc(RANK,...) \
  do { if((RANK)==log_mpid_) log_general_(LOG_STD,  __FILE__, __LINE__, ##__VA_ARGS__); } while (0)

/* Output warnings (prefix WARN) only on given RANK */
#define log_warnc(RANK,...) \
  do { if((RANK)==log_mpid_) log_general_(LOG_WARN,  __FILE__, __LINE__, ##__VA_ARGS__); } while (0)

/*! Initialize logging. Should be called after MPI_Init(), where applicable.
 *  @param fp Either a valid file pointer, or NULL to use stdout/stderr.
 */
void log_init(FILE * fp);

/*! Set output stream.
 *  @param fp Either a valid file pointer, or NULL to use stdout/stderr.
 */
void log_set_output(FILE * fp);

/*! Set logging level.
 *  @param level Minimum level to actually log.
 *  @note Errors and fatal messages are always logged.
 */
void log_set_level(log_Level level);

/*! Set logging level from string.
 *  @param s_level Minimum level to actually log as string.
 *  @note Errors and fatal messages are always logged.
 */
void log_set_level_from_string(const char *s_level);

/*! Get the current logging level.
 *  @return Currently stored logging level.
 */
log_Level log_get_level();

/*! Get a logging level as human-readable character string.
 *  @param level Logging level to query.
 *  @return Logging level as character string.
 */
const char *log_get_level_string(log_Level level);

/*! Get the name of the node an (MPI) process is running on.
 *  @return Character string containing node identifier (usually its name).
 */
const char *log_get_nodename();

/*! Get the identifier of the current (MPI) process.
 *  @return MPI identifier as integer.
 */
int log_get_mpid();

/* 
 * Internal function, should not be called directly. Use logging macros given above.
 */
void log_general_(log_Level level, const char *file, int line, const char *fmt, ...);

/*! Print a banner. Should be called after log_init().
 *  @param prog Main program for which logging is used.
 */
void log_banner(log_Program prog);

/*! Finalize logging and print timing statistics. */
void log_finalize();

#endif                          // __LOGGING_H__
