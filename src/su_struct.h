
#ifndef __SU_STRUCT_H__
#define __SU_STRUCT_H__

typedef struct {
    int tracl;                  /*!< trace sequence number within line (1-4) */
    int tracr;                  /*!< trace sequence number within reel (5-8) */
    int fldr;                   /*!< field record number (9-12) */
    int tracf;                  /*!< trace number within field record (13-16) */
    int ep;                     /*!< energy source point number (17-20) */
    int cdp;                    /*!< ensemble number (21-24) */
    int cdpt;                   /*!< trace number within ensemble (25-28) */
    short trid;                 /*!< trace identification code (29-30) */
    short nvs;                  /*!< number of vertically summed traces (31-32) */
    short nhs;                  /*!< number of horizontally summed traces (33-34) */
    short duse;                 /*!< data use: 1 = production; 2 = test (35-36) */
    int offset;                 /*!< distance from source point to receiver group (37-40) */
    int gelev;                  /*!< receiver group elevation from sea level (41-44) */
    int selev;                  /*!< source elevation from sea level (45-48) */
    int sdepth;                 /*!< source depth (positive) (49-52) */
    int gdel;                   /*!< datum elevation at receiver group (53-56) */
    int sdel;                   /*!< datum elevation at source (57-60) */
    int swdep;                  /*!< water depth at source (61-64) */
    int gwdep;                  /*!< water depth at receiver group (65-68) */
    short scalel;               /*!< scale factor for previous 7 entries (69-70) */
    short scalco;               /*!< scale factor for next 4 entries (71-72) */
    int sx;                     /*!< X source coordinate (73-76) */
    int sy;                     /*!< Y source coordinate (77-80) */
    int gx;                     /*!< X group coordinate (81-84) */
    int gy;                     /*!< Y group coordinate (85-88) */
    short counit;               /*!< coordinate units code (89-90) */
    short wevel;                /*!< weathering velocity (91-92) */
    short swevel;               /*!< subweathering velocity (93-94) */
    short sut;                  /*!< uphole time at source (95-96) */
    short gut;                  /*!< uphole time at receiver group (97-98) */
    short sstat;                /*!< source static correction (99-100) */
    short gstat;                /*!< group static correction (101-102) */
    short tstat;                /*!< total static applied (103-104) */
    short laga;                 /*!< lag time A (105-106) */
    short lagb;                 /*!< lag time B (107-108) */
    short delrt;                /*!< delay recording time (milli-seconds) (109-110) */
    short muts;                 /*!< mute time start (milli-seconds) (111-112) */
    short mute;                 /*!< mute time end (milli-seconds) (113-114) */
    unsigned short ns;          /*!< number of samples in this trace (115-116) */
    unsigned short dt;          /*!< sample interval (micro-seconds) (117-118) */
    short gain;                 /*!< gain type of field instruments code (119-120) */
    short igc;                  /*!< instrument gain constant (dB) (121-122) */
    short igi;                  /*!< instrument early or initial gain (db) (123-124) */
    short corr;                 /*!< correlated: 1 = no ; 2 = yes (125-126) */
    short sfs;                  /*!< sweep frequency at start (Hz) (127-128) */
    short sfe;                  /*!< sweep frequency at end (Hz) (129-130) */
    short slen;                 /*!< sweep length (milli-seconds) (131-132)  */
    short styp;                 /*!< sweep type code (133-134) */
    short stas;                 /*!< sweep trace length at start (milli-seconds) (135-136) */
    short stae;                 /*!< sweep trace length at end (milli-seconds) (137-138) */
    short tatyp;                /*!< taper type: 1=linear, 2=cos^2, 3=other (139-140) */
    short afilf;                /*!< alias filter frequency if used (Hz) (141-142) */
    short afils;                /*!< alias filter slope (dB/octave) (143-144) */
    short nofilf;               /*!< notch filter frequency if used (Hz) (145-146) */
    short nofils;               /*!< notch filter slope (dB/octave) (147-148) */
    short lcf;                  /*!< low cut frequency if used (Hz) (149-150) */
    short hcf;                  /*!< high cut frequncy if used (Hz) (151-152) */
    short lcs;                  /*!< low cut slope (dB/octave) (153-154) */
    short hcs;                  /*!< high cut slope (dB/octave) (155-156) */
    short year;                 /*!< year data recorded (4 digits) (157-158) */
    short day;                  /*!< day of year (159-160) */
    short hour;                 /*!< hour of day (24 hour clock) (161-162) */
    short minute;               /*!< minute of hour (163-164) */
    short sec;                  /*!< second of minute (165-166) */
    short timbas;               /*!< time basis code (167-168) */
    short trwf;                 /*!< trace weighting factor (169-170) */
    short grnors;               /*!< geophone number of roll switch pos one (171-172) */
    short grnofr;               /*!< geophone number of trace one within fldr (173-174) */
    short grnlof;               /*!< geophone number of last trace within fldr (175-176) */
    short gaps;                 /*!< gap size (total number of groups dropped) (177-178) */
    short otrav;                /*!< overtravel taper code (179-180) */
    float d1;                   /*!< sample spacing for non-seismic data */
    float f1;                   /*!< first sample location for non-seismic data */
    float d2;                   /*!< sample spacing between traces */
    float f2;                   /*!< first trace location */
    float ungpow;               /*!< neg of power used for dynamic range compression */
    float unscale;              /*!< reciprocal of scaling factor to normalize range */
    int ntr;                    /*!< number of traces */
    short mark;                 /*!< mark selected traces */
    short unass[15];
} SUhead;

/*! Initialize a SU header structure with zeros.
 *  @param[in] header SU trace header
 */
void init_SUhead(SUhead * header);

#endif
