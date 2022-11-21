/*------------------------------------------------------------------------
 * Copyright (C) 2011 For the list of authors, see file AUTHORS.
 *
 * This file is part of SOFI2D.
 * 
 * SOFI2D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 * 
 * SOFI2D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with SOFI2D. See file COPYING and/or 
  * <http://www.gnu.org/licenses/gpl-2.0.html>.
--------------------------------------------------------------------------*/
/*------------------------------------------------------------------------
 *   program SOFI2D, reading input-parameters from input-file or stdin
 *   that are formatted according to the json standard
 *
 *  ----------------------------------------------------------------------*/

#include <unistd.h>
#include "fd.h"

void read_par_json(char *fileinp, GlobVar *gv) 
{
  /* declaration of extern variables */
  //extern float *FL;
  
  /* definition of local variables */
  int number_readobjects=0,fserr=0;
  char errormessage[STRING_SIZE];

  //allocate first object in list
  char **varname_list = malloc(STRING_SIZE*sizeof(char*));
  char **value_list = malloc(STRING_SIZE*sizeof(char*));

  //read in objects from file
  number_readobjects=read_objects_from_intputfile(gv->FP, fileinp, varname_list, value_list);
  fprintf(gv->FP,"\nFrom input file %s, %i objects have been read in. \n",fileinp, number_readobjects);

  //print objects to screen
  fprintf(gv->FP, "\n===========================================================");
  fprintf(gv->FP, "\n=   List of Parameters read by the built in Json Parser   =");
  print_objectlist_screen(gv->FP, number_readobjects, varname_list, value_list);

  //extract variables form object list

  /*=================================
    section general grid and discretization parameters
    =================================*/

  if (get_int_from_objectlist("NPROCX",number_readobjects,&(gv->NPROCX),varname_list, value_list))
    declare_error("Variable NPROCX could not be retrieved from the json input file!");
  if (get_int_from_objectlist("NPROCY",number_readobjects,&(gv->NPROCY),varname_list, value_list))
    declare_error("Variable NPROCY could not be retrieved from the json input file!");
  gv->NPROC = (gv->NPROCX)*(gv->NPROCY);
  /*	if (get_int_from_objectlist("RSG",number_readobjects,&(gv->RSG),varname_list, value_list))
	declare_error("Variable RSG could not be retrieved from the json input file!");*/
  if (get_int_from_objectlist("FDORDER",number_readobjects,&(gv->FDORDER),varname_list, value_list))
    declare_error("Variable FDORDER could not be retrieved from the json input file!");
  if (get_int_from_objectlist("FDORDER_TIME",number_readobjects,&(gv->FDORDER_TIME),varname_list, value_list)) {
    (gv->FDORDER_TIME)=2;
  } else {
    if((gv->FDORDER_TIME)!=2 && (gv->FDORDER_TIME)!=4) {
      declare_error("Only FDORDER_TIME 2 or 4 are supported!");
    }
  }
  if (get_int_from_objectlist("MAXRELERROR",number_readobjects,&(gv->MAXRELERROR),varname_list, value_list))
    declare_error("Variable MAXRELERROR could not be retrieved from the json input file!");
  if (get_int_from_objectlist("NX",number_readobjects,&(gv->NX),varname_list, value_list))
    declare_error("Variable NX could not be retrieved from the json input file!");
  if (get_int_from_objectlist("NY",number_readobjects,&(gv->NY),varname_list, value_list))
    declare_error("Variable NY could not be retrieved from the json input file!");
  if (get_float_from_objectlist("DH",number_readobjects,&(gv->DH),varname_list, value_list))
    declare_error("Variable DH could not be retrieved from the json input file!");
  if (get_float_from_objectlist("TIME",number_readobjects,&(gv->TIME),varname_list, value_list))
    declare_error("Variable TIME could not be retrieved from the json input file!");
  if (get_float_from_objectlist("DT",number_readobjects,&(gv->DT),varname_list, value_list))
    declare_error("Variable DT could not be retrieved from the json input file!");
  
  if (get_int_from_objectlist("WEQ",number_readobjects,&(gv->WEQ),varname_list, value_list)) {
    (gv->WEQ)=3;
  } else {
    if((gv->WEQ)!=3 && (gv->WEQ)!=4 && (gv->WEQ)!=5 && (gv->WEQ)!=6 && (gv->WEQ)!=7 && (gv->WEQ)!=8){
      fprintf(gv->FP, " Only the following wave equations are supported: \n");
      fprintf(gv->FP, " WEQ=3 (isotropic elastic) \n");
      fprintf(gv->FP, " WEQ=4 (isotropic viscoelastic) \n");
      fprintf(gv->FP, " WEQ=5 (elastic VTI) \n");
      fprintf(gv->FP, " WEQ=6 (viscoelastic VTI) \n");
      fprintf(gv->FP, " WEQ=7 (elastic TTI) \n");
      fprintf(gv->FP, " WEQ=8 (viscoelastic TTI) \n");
      declare_error("Stop.");
    }
  }
  
  /*=================================
    section source parameters
    =================================*/
  
  if (get_int_from_objectlist("SOURCE_TYPE",number_readobjects,&(gv->SOURCE_TYPE),varname_list, value_list))
    declare_error("Variable SOURCE_TYPE could not be retrieved from the json input file!");
  if (get_int_from_objectlist("SIGOUT",number_readobjects,&(gv->SIGOUT),varname_list, value_list)) {
      //declare_error("Variable SIGOUT could not be retrieved from the json input file!");
      fprintf(gv->FP, " ADDITIONAL VARIABLE OPTIONS: With the variables SIGOUT, SIGOUT_FILE & SIGOUT_FORMAT you can output the source signal to file. \n\n\n");
  } else {
    if ((gv->SIGOUT)==1) {
      if (get_string_from_objectlist("SIGOUT_FILE",number_readobjects,gv->SIGOUT_FILE,varname_list, value_list)) {
          declare_error("Variable SIGOUT_FILE could not be retrieved from the json input file!");
      }
      if (get_int_from_objectlist("SIGOUT_FORMAT",number_readobjects,&(gv->SIGOUT_FORMAT),varname_list, value_list)) {
        declare_error("Variable SIGOUT_FORMAT could not be retrieved from the json input file!");
      }
    }
  }
  if (get_int_from_objectlist("SOURCE_SHAPE",number_readobjects,&(gv->SOURCE_SHAPE),varname_list, value_list))
    declare_error("Variable SOURCE_SHAPE could not be retrieved from the json input file!");
  else {
    if ((gv->SOURCE_SHAPE)==3) {
      if (get_string_from_objectlist("SIGNAL_FILE",number_readobjects,(gv->SIGNAL_FILE),varname_list, value_list))
	declare_error("Variable SIGNAL_FILE could not be retrieved from the json input file!");
    }
  }
  if (get_int_from_objectlist("SRCREC",number_readobjects,&(gv->SRCREC),varname_list, value_list))
    declare_error("Variable SRCREC could not be retrieved from the json input file!");
  else {
    if ((gv->SRCREC)==1) {
      if (get_string_from_objectlist("SOURCE_FILE",number_readobjects,(gv->SOURCE_FILE),varname_list, value_list))
	declare_error("Variable SOURCE_FILE could not be retrieved from the json input file!");
      if (get_int_from_objectlist("RUN_MULTIPLE_SHOTS",number_readobjects,&(gv->RUN_MULTIPLE_SHOTS),varname_list, value_list))
	declare_error("Variable RUN_MULTIPLE_SHOTS could not be retrieved from the json input file!");
    }
    if ((gv->SRCREC)==2) {
      if (get_float_from_objectlist("PLANE_WAVE_DEPTH",number_readobjects,&(gv->PLANE_WAVE_DEPTH),varname_list, value_list))
	declare_error("Variable PLANE_WAVE_DEPTH could not be retrieved from the json input file!");
      else {
	if ((gv->PLANE_WAVE_DEPTH)>0.0) {
	  if (get_float_from_objectlist("PLANE_WAVE_ANGLE",number_readobjects,&(gv->PLANE_WAVE_ANGLE),varname_list, value_list))
	    declare_error("Variable PLANE_WAVE_ANGLE could not be retrieved from the json input file!");
	  if (get_float_from_objectlist("TS",number_readobjects,&(gv->TS),varname_list, value_list))
	    declare_error("Variable TS could not be retrieved from the json input file!");
	}
      }
    }
  }
  
  /*=================================
    section boundary parameters
    =================================*/
  
  if (get_int_from_objectlist("FREE_SURF",number_readobjects,&(gv->FREE_SURF),varname_list, value_list))
    declare_error("Variable FREE_SURF could not be retrieved from the json input file!");
  if (get_int_from_objectlist("BOUNDARY",number_readobjects,&(gv->BOUNDARY),varname_list, value_list))
    declare_error("Variable BOUNDARY could not be retrieved from the json input file!");
  if (get_int_from_objectlist("FW",number_readobjects,&(gv->FW),varname_list, value_list))
    declare_error("Variable FW could not be retrieved from the json input file!");
  if (get_int_from_objectlist("ABS_TYPE",number_readobjects,&(gv->ABS_TYPE),varname_list, value_list))
    declare_error("Variable ABS_TYPE could not be retrieved from the json input file!");
  
  if ((gv->ABS_TYPE)==1) {
    if (get_float_from_objectlist("NPOWER",number_readobjects,&(gv->NPOWER),varname_list, value_list)) 
      declare_error("Variable NPOWER could not be retrieved from the json input file!");
    if (get_float_from_objectlist("K_MAX_CPML",number_readobjects,&(gv->K_MAX_CPML),varname_list, value_list)) 
      declare_error("Variable K_MAX_CPML could not be retrieved from the json input file!");
    if (get_float_from_objectlist("FPML",number_readobjects,&(gv->FPML),varname_list, value_list))
      declare_error("Variable FPML could not be retrieved from the json input file!");
    if (get_float_from_objectlist("VPPML",number_readobjects,&(gv->VPPML),varname_list, value_list))
      declare_error("Variable VPPML could not be retrieved from the json input file!");
  }
  if ((gv->ABS_TYPE)==2) {
    if (get_float_from_objectlist("DAMPING",number_readobjects,&(gv->DAMPING),varname_list, value_list))
      declare_error("Variable DAMPING could not be retrieved from the json input file!");
  }
  
  /*=================================
    section snapshot parameters
    =================================*/

  if (get_int_from_objectlist("SNAP",number_readobjects,&(gv->SNAP),varname_list, value_list))
    declare_error("Variable SNAP could not be retrieved from the json input file!");
  else {
    if ((gv->SNAP)>0) {
      if (get_int_from_objectlist("SNAP_FORMAT",number_readobjects,&(gv->SNAP_FORMAT),varname_list, value_list))
	declare_error("Variable SNAP_FORMAT could not be retrieved from the json input file!");
      if (get_float_from_objectlist("TSNAP1",number_readobjects,&(gv->TSNAP1),varname_list, value_list))
	declare_error("Variable TSNAP1 could not be retrieved from the json input file!");
      if (get_float_from_objectlist("TSNAP2",number_readobjects,&(gv->TSNAP2),varname_list, value_list))
	declare_error("Variable TSNAP2 could not be retrieved from the json input file!");
      if (get_float_from_objectlist("TSNAPINC",number_readobjects,&(gv->TSNAPINC),varname_list, value_list))
	declare_error("Variable TSNAPINC could not be retrieved from the json input file!");
      if (get_string_from_objectlist("SNAP_FILE",number_readobjects,(gv->SNAP_FILE),varname_list, value_list))
	declare_error("Variable SNAP_FILE could not be retrieved from the json input file!");
    }
  }
  /* increments are read in any case, because they will be also used as increment for model output */
  if (get_int_from_objectlist("IDX",number_readobjects,&(gv->IDX),varname_list, value_list))
    declare_error("Variable IDX could not be retrieved from the json input file!");
  if (get_int_from_objectlist("IDY",number_readobjects,&(gv->IDY),varname_list, value_list))
    declare_error("Variable IDY could not be retrieved from the json input file!");
  
  /*=================================
    section seismogramm parameters
    =================================*/
  
  if (get_int_from_objectlist("SEISMO",number_readobjects,&(gv->SEISMO),varname_list, value_list))
    declare_error("Variable SEISMO could not be retrieved from the json input file!");
  else {
    if ((gv->SEISMO)>0){
      if (get_string_from_objectlist("REC_FILE",number_readobjects,(gv->REC_FILE),varname_list, value_list))
	declare_error("Variable REC_FILE could not be retrieved from the json input file!");
      if (get_string_from_objectlist("SEIS_FILE",number_readobjects,(gv->SEIS_FILE),varname_list, value_list))
	declare_error("Variable SEIS_FILE could not be retrieved from the json input file!");
      if (get_int_from_objectlist("READREC",number_readobjects,&(gv->READREC),varname_list, value_list))
	declare_error("Variable READREC could not be retrieved from the json input file!");
      else {
	if ((gv->READREC)==0) {
	  if (get_float_from_objectlist("XREC1",number_readobjects,&(gv->XREC1),varname_list, value_list))
	    declare_error("Variable XREC1 could not be retrieved from the json input file!");
	  if (get_float_from_objectlist("XREC2",number_readobjects,&(gv->XREC2),varname_list, value_list))
	    declare_error("Variable XREC2T could not be retrieved from the json input file!");
	  if (get_float_from_objectlist("YREC1",number_readobjects,&(gv->YREC1),varname_list, value_list))
	    declare_error("Variable YREC1 could not be retrieved from the json input file!");
	  if (get_float_from_objectlist("YREC2",number_readobjects,&(gv->YREC2),varname_list, value_list))
	    declare_error("Variable YREC2 could not be retrieved from the json input file!");
	}
      }
      if (get_int_from_objectlist("NDT",number_readobjects,&(gv->NDT),varname_list, value_list))
	declare_error("Variable NDT could not be retrieved from the json input file!");
      if (get_int_from_objectlist("SEIS_FORMAT",number_readobjects,&(gv->SEIS_FORMAT),varname_list, value_list))
	declare_error("Variable SEIS_FORMAT could not be retrieved from the json input file!");
      
      if (get_int_from_objectlist("REC_ARRAY",number_readobjects,&(gv->REC_ARRAY),varname_list, value_list))
	declare_error("Variable REC_ARRAY could not be retrieved from the json input file!");
      else {
	if ((gv->REC_ARRAY)>0) {
	  if (get_int_from_objectlist("DRX",number_readobjects,&(gv->DRX),varname_list, value_list))
	    declare_error("Variable DRX could not be retrieved from the json input file!");
	  if (get_float_from_objectlist("REC_ARRAY_DEPTH",number_readobjects,&(gv->REC_ARRAY_DEPTH),varname_list, value_list))
	    declare_error("Variable REC_ARRAY_DEPTH could not be retrieved from the json input file!");
	  if (get_float_from_objectlist("REC_ARRAY_DIST",number_readobjects,&(gv->REC_ARRAY_DIST),varname_list, value_list))
	    declare_error("Variable REC_ARRAY_DIST could not be retrieved from the json input file!");
	}
      }
      if (get_float_from_objectlist("REFRECX",number_readobjects,&(gv->REFREC[1]),varname_list, value_list))
	declare_error("Variable REFRECX could not be retrieved from the json input file!");
      if (get_float_from_objectlist("REFRECY",number_readobjects,&(gv->REFREC[2]),varname_list, value_list))
	declare_error("Variable REFRECY could not be retrieved from the json input file!");
      if (get_float_from_objectlist("NGEOPH",number_readobjects,&(gv->NGEOPH),varname_list, value_list))
	declare_error("Variable NGEOPH could not be retrieved from the json input file!");
    }
  }

  /*=================================
    section general model and log parameters
    =================================*/
  if (get_string_from_objectlist("MFILE",number_readobjects,(gv->MFILE),varname_list, value_list))
    declare_error("Variable MFILE could not be retrieved from the json input file!");
  if (get_int_from_objectlist("WRITE_MODELFILES",number_readobjects,&(gv->WRITE_MODELFILES),varname_list, value_list)) {}
  if (get_int_from_objectlist("LOG",number_readobjects,&(gv->LOG),varname_list, value_list))
    declare_error("Variable LOG could not be retrieved from the json input file!");
  if (get_int_from_objectlist("CHECKPTREAD",number_readobjects,&(gv->CHECKPTREAD),varname_list, value_list))
    declare_error("Variable CHECKPTREAD could not be retrieved from the json input file!");
  if (get_int_from_objectlist("CHECKPTWRITE",number_readobjects,&(gv->CHECKPTWRITE),varname_list, value_list))
    declare_error("Variable CHECKPTWRITE could not be retrieved from the json input file!");
  if (get_string_from_objectlist("LOG_FILE",number_readobjects,(gv->LOG_FILE),varname_list, value_list))
    declare_error("Variable LOG_FILE could not be retrieved from the json input file!");
  if (get_string_from_objectlist("CHECKPT_FILE",number_readobjects,(gv->CHECKPTFILE),varname_list, value_list))
    declare_error("Variable CHECKPT_FILE could not be retrieved from the json input file!");
  if (get_int_from_objectlist("READMOD",number_readobjects,&(gv->READMOD),varname_list, value_list))
    declare_error("Variable READMOD could not be retrieved from the json input file!");
  if (get_int_from_objectlist("OUT_TIMESTEP_INFO",number_readobjects,&(gv->OUTNTIMESTEPINFO),varname_list, value_list)) {}
  
  if (get_float_from_objectlist("TAU",number_readobjects,&(gv->TAU),varname_list, value_list))
    declare_error("Variable TAU could not be retrieved from the json input file!");
  if (get_int_from_objectlist("L",number_readobjects,&(gv->L),varname_list, value_list))
    declare_error("Variable L could not be retrieved from the json input file!");
  
  if (((gv->L)==0) && (((gv->WEQ)==4) || ((gv->WEQ)==6) || ((gv->WEQ)==8)))
    declare_error("L>0 for viscoelastic simulations (WEQ=4 or WEQ=6 or WEQ=8)!");
  
  if (((gv->L)>0) && (((gv->WEQ)==3) || ((gv->WEQ)==5) || ((gv->WEQ)==7))) {
    fprintf(gv->FP, "L is set to zero for elastic simulation (WEQ=3 or WEQ=5 or WEQ=7) ! \n\n\n");
    (gv->L) = 0;
  }
  
  if ((gv->L)>0) {
    (gv->FL) = vector(1,(gv->L));
    switch((gv->L)) {
    case 0:
      break;
    case 5:
      if (get_float_from_objectlist("FL5",number_readobjects,&(gv->FL[5]),varname_list, value_list))
	declare_error("Variable FL5 could not be retrieved from the json input file!");
    case 4:
      if (get_float_from_objectlist("FL4",number_readobjects,&(gv->FL[4]),varname_list, value_list))
	declare_error("Variable FL4 could not be retrieved from the json input file!");
    case 3:
      if (get_float_from_objectlist("FL3",number_readobjects,&(gv->FL[3]),varname_list, value_list))
	declare_error("Variable FL3 could not be retrieved from the json input file!");
    case 2:
      if (get_float_from_objectlist("FL2",number_readobjects,&(gv->FL[2]),varname_list, value_list))
	declare_error("Variable FL2 could not be retrieved from the json input file!");
    case 1:
      if (get_float_from_objectlist("FL1",number_readobjects,&(gv->FL[1]),varname_list, value_list))
	declare_error("Variable FL1 could not be retrieved from the json input file!");
      break;
    default:
      declare_error("More than five relaxation parameters (L>5) not implemented yet!");
      break;
    }
  } else {
    (gv->FL) = NULL;
  }
  
  /********************************************/
  /* Check files and directories if necessary */
  /********************************************/
  
  /* signal file */
  if ((gv->SOURCE_SHAPE) == 3) {
    if (access((gv->SIGNAL_FILE),0) != 0) {
      fprintf(gv->FP, "\n==================================================================\n");
      fprintf(gv->FP, "  ERROR parsing input file <%s>:\n", fileinp);
      fprintf(gv->FP, "        The signal file does not exist!\n");
      fprintf(gv->FP, "        File name: <%s>", (gv->SIGNAL_FILE));
      fprintf(gv->FP, "\n==================================================================\n");
      fserr = 1;
    } else if (access((gv->SIGNAL_FILE),4) != 0) {
      fprintf(gv->FP, "\n==================================================================\n");
      fprintf(gv->FP, "  ERROR parsing input file <%s>:\n", fileinp);
      fprintf(gv->FP, "        The signal file does not have read access!\n");
      fprintf(gv->FP, "        File name: <%s>", (gv->SIGNAL_FILE));
      fprintf(gv->FP, "\n==================================================================\n");
      fserr = 1;
    }
  }
  
  /* source file */
  if ((gv->SRCREC)==1) {
    if (access((gv->SOURCE_FILE),0) != 0) {
      fprintf(gv->FP, "\n==================================================================\n");
      fprintf(gv->FP, "  ERROR parsing input file <%s>:\n", fileinp);
      fprintf(gv->FP, "        The source file does not exist!\n");
      fprintf(gv->FP, "        File name: <%s>", (gv->SOURCE_FILE));
      fprintf(gv->FP, "\n==================================================================\n");
      fserr = 1;
    } else if (access((gv->SOURCE_FILE),4) != 0) {
      fprintf(gv->FP, "\n==================================================================\n");
      fprintf(gv->FP, "  ERROR parsing input file <%s>:\n", fileinp);
      fprintf(gv->FP, "        The source file does not have read access!\n");
      fprintf(gv->FP, "        File name: <%s>", (gv->SOURCE_FILE));
      fprintf(gv->FP, "\n==================================================================\n");
      fserr = 1;
    }
  }
  
  /* receiver file */
  if ((gv->READREC)) {
    if (access((gv->REC_FILE),0) != 0) {
      fprintf(gv->FP, "\n==================================================================\n");
      fprintf(gv->FP, "  ERROR parsing input file <%s>:\n", fileinp);
      fprintf(gv->FP, "        The receiver file does not exist!\n");
      fprintf(gv->FP, "        File name: <%s>", (gv->REC_FILE));
      fprintf(gv->FP, "\n==================================================================\n");
      fserr = 1;
    } else if (access((gv->REC_FILE),4) != 0) {
      fprintf(gv->FP, "\n==================================================================\n");
      fprintf(gv->FP, "  ERROR parsing input file <%s>:\n", fileinp);
      fprintf(gv->FP, "        The receiver file does not have read access!\n");
      fprintf(gv->FP, "        File name: <%s>", (gv->REC_FILE));
      fprintf(gv->FP, "\n==================================================================\n");
      fserr = 1;
    }
  }
  
  /* checkpoint file */
  if ((gv->CHECKPTREAD) || (gv->CHECKPTWRITE)) {
    if (access((gv->CHECKPTFILE),0) != 0) {
      fprintf(gv->FP, "\n==================================================================\n");
      fprintf(gv->FP, "  ERROR parsing input file <%s>:\n", fileinp);
      fprintf(gv->FP, "        The checkpoint file does not exist!\n");
      fprintf(gv->FP, "        File name: <%s>", (gv->CHECKPTFILE));
      fprintf(gv->FP, "\n==================================================================\n");
      fserr = 1;
    } else if (access((gv->CHECKPTFILE),6) != 0) {
      fprintf(gv->FP, "\n==================================================================\n");
      fprintf(gv->FP, "  ERROR parsing input file <%s>:\n", fileinp);
      fprintf(gv->FP, "        The checkpoint file does not have read and/or write access!\n");
      fprintf(gv->FP, "        File name: <%s>", (gv->CHECKPTFILE));
      fprintf(gv->FP, "\n==================================================================\n");
      fserr = 1;
    }
  }
  
  /********************************************/
  /* ERROR                                    */
  /********************************************/
  if (fserr) {
    fprintf(gv->FP, "\n");
    sprintf(errormessage, "\n  in: <read_par_json.c> \n");
    declare_error(errormessage);
  }
}
