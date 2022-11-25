/*------------------------------------------------------------------------
 * Copyright (C) 2016 For the list of authors, see file AUTHORS.
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
/* -------------------------------------------------------------
 * function json_parser
 * reads input parameters from input file formatted according to the json standard
 * also see www.json.org
 *
 *------------------------------------------------------------- */

#include "fd.h"
#include "logging.h"

int read_objects_from_intputfile(const char *input_file,char **varname_list,char **value_list) 
{
  char varname_tmp1[STRING_SIZE], varname_tmp2[STRING_SIZE], varname_tmp3[STRING_SIZE];
  char varname_tmp4[STRING_SIZE], varname_tmp5[STRING_SIZE];
  char value_tmp1[STRING_SIZE], value_tmp2[STRING_SIZE], value_tmp3[STRING_SIZE];
  char value_tmp4[STRING_SIZE], value_tmp5[STRING_SIZE];
  char cline[STRING_SIZE];
  int occurence_doublequotes = 0, occurence_commas = 0;
  int lineno=0;
  int number_readobject=0;
  FILE *fp_in = NULL;
  
  //Open parameter input file
  fp_in=fopen(input_file,"r");
  if (!fp_in) {
    log_fatal("Could not open file %s for reading.\n", input_file);
  }
  
  //read line by line into a string covering the whole line
  while (fgets(cline,STRING_SIZE,fp_in)) { /* leaves via break */
    /* If there are more than 255 characters in one line, this does not work. */
    //count of line numbers
    lineno++;
    
    /* tests if line is NOT a comment line*/
    /* tests if line contains at least a colon, double quote per line*/
    if (((strstr(cline,":"))&&((strstr(cline,"\"")))) && (!(strstr(cline,"comment")) && !(strstr(cline,"Comment"))))  {
      
      //count number of double quotes and colon signs
      occurence_doublequotes=count_occure_charinstring(cline,"\"");
      occurence_commas=count_occure_charinstring(cline,",");
      
      //only two pais of double quotes are allowed per line
      switch (occurence_doublequotes) {
      case 4:
	//up to 5 objects can be defined per line, more can be implemented here
	switch (occurence_commas) {
	case 0: //only a single object (name+value) in line
		//remove old data from strings
	  memset(value_tmp1, '\0', sizeof(value_tmp1));
	  memset(varname_tmp1, '\0', sizeof(varname_tmp1));
	  
	  //extract object name + object value from the line-string
	  if (sscanf(cline," \"%[^\"]\" : \"%[^\"]\"",varname_tmp1,value_tmp1) != 2) {
	    log_fatal("Error in input file, line %i, cannot read object name and object value!\n",lineno);
	  }
	  //add extracted strings to object list
	  add_object_tolist(varname_tmp1, value_tmp1,&number_readobject, varname_list, value_list);
	  break;
	case 1: //only a single object (name+value) in line
		//remove old data from strings
	  memset(value_tmp1, '\0', sizeof(value_tmp1));
	  memset(varname_tmp1, '\0', sizeof(varname_tmp1));
	  //extract object name + object value from the line-string
	  if (sscanf(cline," \"%[^\"]\" : \"%[^\"]\"",varname_tmp1,value_tmp1) != 2) {
	    log_fatal("Error in input file, line %i, cannot read object name and object value!\n",lineno);
	  }
	  //add extracted strings to object list
	  add_object_tolist(varname_tmp1, value_tmp1,&number_readobject, varname_list, value_list);
	  break;
	case 3://two objects (name+value) in line
	  //remove old data from strings
	  memset(value_tmp1, '\0', sizeof(value_tmp1));
	  memset(varname_tmp1, '\0', sizeof(varname_tmp1));
	  memset(value_tmp2, '\0', sizeof(value_tmp2));
	  memset(varname_tmp2, '\0', sizeof(varname_tmp2));
	  //extract object name + object value from the line-string
	  if (sscanf(cline," \"%[^,],%[^\"]\" : \"%[^,],%[^\"]\"",
		     varname_tmp1,varname_tmp2,value_tmp1,value_tmp2) != 4) {
	    log_fatal("Error in input file, line %i, cannot read two object names and values!\n",lineno);
	  }
	  //add extracted strings to object list
	  add_object_tolist(varname_tmp1, value_tmp1,&number_readobject, varname_list, value_list);
	  add_object_tolist(varname_tmp2, value_tmp2,&number_readobject, varname_list, value_list);
	  break;
	case 5://three objects (name+value) in line
	  //remove old data from strings
	  memset(value_tmp1, '\0', sizeof(value_tmp1));
	  memset(value_tmp2, '\0', sizeof(value_tmp2));
	  memset(value_tmp3, '\0', sizeof(value_tmp3));
	  memset(varname_tmp1, '\0', sizeof(varname_tmp1));
	  memset(varname_tmp2, '\0', sizeof(varname_tmp2));
	  memset(varname_tmp3, '\0', sizeof(varname_tmp3));
	  if (sscanf(cline," \"%[^,],%[^,],%[^\"]\" : \"%[^,],%[^,],%[^\"]\"",
		     varname_tmp1,varname_tmp2,varname_tmp3,value_tmp1,value_tmp2,value_tmp3) != 6) {
	    log_fatal("Error in input file, line %i, cannot read three object names and values!\n",lineno);
	  }
	  add_object_tolist(varname_tmp1, value_tmp1,&number_readobject, varname_list, value_list);
	  add_object_tolist(varname_tmp2, value_tmp2,&number_readobject, varname_list, value_list);
	  add_object_tolist(varname_tmp3, value_tmp3,&number_readobject, varname_list, value_list);
	  break;
	case 7://four objects (name+value) in line
	  //remove old data from strings
	  memset(value_tmp1, '\0', sizeof(value_tmp1));
	  memset(value_tmp2, '\0', sizeof(value_tmp2));
	  memset(value_tmp3, '\0', sizeof(value_tmp3));
	  memset(value_tmp4, '\0', sizeof(value_tmp4));
	  memset(varname_tmp1, '\0', sizeof(varname_tmp1));
	  memset(varname_tmp2, '\0', sizeof(varname_tmp2));
	  memset(varname_tmp3, '\0', sizeof(varname_tmp3));
	  memset(varname_tmp4, '\0', sizeof(varname_tmp4));
	  if (sscanf(cline," \"%[^,],%[^,],%[^,],%[^\"]\" : \"%[^,],%[^,],%[^,],%[^\"]\"",
		     varname_tmp1,varname_tmp2,varname_tmp3,varname_tmp4,
		     value_tmp1,value_tmp2,value_tmp3,value_tmp4) != 8) {
	    log_fatal("Error in input file, line %i, cannot read three object names and values!\n",lineno);
	  }
	  add_object_tolist(varname_tmp1, value_tmp1,&number_readobject, varname_list, value_list);
	  add_object_tolist(varname_tmp2, value_tmp2,&number_readobject, varname_list, value_list);
	  add_object_tolist(varname_tmp3, value_tmp3,&number_readobject, varname_list, value_list);
	  add_object_tolist(varname_tmp4, value_tmp4,&number_readobject, varname_list, value_list);
	  break;
	case 9://five objects (name+value) in line
	  //remove old data from strings
	  memset(value_tmp1, '\0', sizeof(value_tmp1));
	  memset(value_tmp2, '\0', sizeof(value_tmp2));
	  memset(value_tmp3, '\0', sizeof(value_tmp3));
	  memset(value_tmp4, '\0', sizeof(value_tmp4));
	  memset(value_tmp5, '\0', sizeof(value_tmp5));
	  memset(varname_tmp1, '\0', sizeof(varname_tmp1));
	  memset(varname_tmp2, '\0', sizeof(varname_tmp2));
	  memset(varname_tmp3, '\0', sizeof(varname_tmp3));
	  memset(varname_tmp4, '\0', sizeof(varname_tmp4));
	  memset(varname_tmp5, '\0', sizeof(varname_tmp5));
	  if (sscanf(cline," \"%[^,],%[^,],%[^,],%[^,],%[^\"]\" : \"%[^,],%[^,],%[^,],%[^,],%[^\"]\"",
		     varname_tmp1,varname_tmp2,varname_tmp3,varname_tmp4,varname_tmp5,
		     value_tmp1,value_tmp2,value_tmp3,value_tmp4,value_tmp5) != 10) {
	    log_fatal("Error in input file, line %i, cannot read three object names and values!\n",lineno);
	  }
	  add_object_tolist(varname_tmp1, value_tmp1,&number_readobject, varname_list, value_list);
	  add_object_tolist(varname_tmp2, value_tmp2,&number_readobject, varname_list, value_list);
	  add_object_tolist(varname_tmp3, value_tmp3,&number_readobject, varname_list, value_list);
	  add_object_tolist(varname_tmp4, value_tmp4,&number_readobject, varname_list, value_list);
	  add_object_tolist(varname_tmp5, value_tmp5,&number_readobject, varname_list, value_list);
	  //very strange: code crashes if both lines are commented here!
	  //this only effects the last case of the switch!
	  //should though not affect anything as long as number_readobject keeps its value
	  //in this case a new object is allocated which is already there...
	  //varname_list = malloc(sizeof(*varname_list));
	  //value_list = malloc(sizeof(*value_list));
	  varname_list[number_readobject] = malloc(STRING_SIZE*sizeof(char *));
	  //varname_list[number_readobject] = (char**)malloc(STRING_SIZE*sizeof(char*));
	  //value_list[number_readobject] = (char**)malloc(STRING_SIZE*sizeof(char*));
	  break;
	default:
	  log_fatal("Error in input file, line %i, only 0, 1, 3, 5, 7 or 9 commas are allowed per line, but found %i!\n",lineno,occurence_commas);
	  break;
	}
	break;
      default:
	log_fatal("Error in input file, line %i, only 4 (two pairs) of double quotes are allowed per line, but found %i!\n",lineno,occurence_doublequotes);
	break;
      }
      //log_std("JSON - line %i contains object # %i varname %s with value %s.\n",
      //        lineno,number_readobject, varname_list[number_readobject-1],value_list[number_readobject-1]);
    }
  }
  fclose(fp_in);
  return number_readobject;
}

void print_objectlist_screen(int number_readobject,char **varname_list,char **value_list) 
{
  /* if we are not on DEBUG level, return straight away, there is no point in
     going through the loop in the first place, it would produce no output */
  if (log_get_level() < LOG_DEBUG) return;
  
  log_debug("==========================================================\n");
  log_debug("   # |        object name | object value                  \n");
  log_debug("==========================================================\n");
  for (int ii=0; ii<number_readobject; ii++) {
    log_debug(" %3.0i | %18s | %s\n",ii+1, varname_list[ii],value_list[ii]);
  }
  log_debug("==========================================================\n");
  
  return;
}

int count_occure_charinstring(char stringline[STRING_SIZE], char teststring[]) 
{
  int number_occurence=0;
  
  for (int ii=0; stringline[ii] != '\0'; ii++) {
    //log_std("lineno = %i ii= %i cline[ii]= %c  ",lineno, ii,cline[ii]);
    //log_std("teststring= %s ",teststring);
    if (strchr(teststring,stringline[ii])) {
      number_occurence++;
    }
  }
  
  return number_occurence;
}

void copy_str2str_uptochar(char string_in[STRING_SIZE], char string_out[STRING_SIZE], char teststring[]) 
{
  for (int ii=0; string_in[ii] != '\0'; ii++) {
    if (strchr(teststring,string_in[ii])) {
      strncpy(string_out,string_in,ii);
    }
  }
  
  return;
}

int get_int_from_objectlist(char string_in[STRING_SIZE], int number_readobject, int *int_buffer,
                            char **varname_list,char **value_list, int *used_list)
{
  int ii=0, checkifstringfound=1;
  double double_buffer;
  double intpart;
  char *string_buffer;
  
  while ((strcmp(varname_list[ii],string_in)!=0) && ((ii+1)<number_readobject)) {
    ii++;
  }
  
  if (strcmp(varname_list[ii],string_in)==0) {
    if (strlen(value_list[ii])==0) {
      log_fatal("Error in input file, value of object %s is empty!\n",string_in);
    }
    memset(&string_buffer, '\0', sizeof(*string_buffer));
    double_buffer = strtod(value_list[ii],&string_buffer);
    if (strlen(string_buffer)>0) {
      /* string empty or 'garbage after double' */
      log_fatal("Error in input file, value of object %s contains more than one float: '%s'!\n",string_in,string_buffer);
    }
    if ((modf(double_buffer, &intpart))==0) {
      *int_buffer = atoi(value_list[ii]);
      checkifstringfound=0; 
    } else {
      //double read, not an int (there are decimal places)
      log_fatal("Error in input file, value of object %s is not an int : %f!\n",string_in,double_buffer);
      //*int_buffer=-1;
      //checkifstringfound=2;
    }  
  } else {
    checkifstringfound=1;
  }

  used_list[ii] = 1;
  
  return checkifstringfound;
}

int get_float_from_objectlist(char string_in[STRING_SIZE], int number_readobject, float *double_buffer,
                              char **varname_list,char **value_list, int *used_list) 
{
  int ii=0, checkifstringfound=1;
  double double_dummy;
  char *string_buffer;
  
  while ((strcmp(varname_list[ii],string_in)!=0) && ((ii+1)<number_readobject)) {
    ii++;
  }
  
  //note: strstr compares if string_in is within varname_list[ii]
  // strcmp compares characterwise if string_in is equal to varname_list[ii]
  if (strcmp(varname_list[ii],string_in)==0) {
    if (strlen(value_list[ii])==0) {
      log_fatal("Error in input file, value of object %s is empty!\n",string_in);
    }
    memset(&string_buffer, '\0', sizeof(*string_buffer));
    double_dummy = strtod(value_list[ii],&string_buffer);
    if ((strlen(string_buffer)==0) || ((strlen(string_buffer)>0) && ((is_string_blankspace(string_buffer))==1))) {
      *double_buffer=double_dummy;
      checkifstringfound=0;
    } else {
      /* string empty or 'garbage after double' */
      log_fatal("Error in input file, value of object %s contains more than one float: '%s'!\n",string_in,string_buffer);
      //checkifstringfound=2;
    }  
  } else {
    checkifstringfound=1;
  }

  used_list[ii] = 1;
  
  return checkifstringfound;
}

int get_string_from_objectlist(char string_in[STRING_SIZE], int number_readobject, char string_buffer[STRING_SIZE],
                               char **varname_list,char **value_list, int *used_list) 
{
  int ii=0, checkifstringfound=1;
  
  while ((strcmp(varname_list[ii],string_in)!=0) && ((ii+1)<number_readobject)) {
    ii++;
  }
  
  if (strcmp(varname_list[ii],string_in)==0) {
    if (strlen(value_list[ii])==0) {
      log_fatal("Error in input file, value of object %s is empty!\n",string_in);
    } else {
      memset(string_buffer, '\0', sizeof(&string_buffer));
      strcpy(string_buffer,value_list[ii]);
      checkifstringfound=0;
    }
  } else {
    checkifstringfound=1;
  }

  used_list[ii] = 1;

  return checkifstringfound;
}

int is_string_blankspace(char string_in[STRING_SIZE]) {

  int ii=0, blank_num=0, string_length=0;

  string_length=strlen(string_in);

  while (ii < string_length)  {
    if (string_in[ii++] == ' ')
      blank_num++;
  }
  
  //log_std("String has length %i and contains %i black spaces.\n",string_length, blank_num);
  if (blank_num==string_length) return 1;
  else return 0;
}

void remove_blankspaces_around_string(char string_in[STRING_SIZE]) 
{
  char string_dummy[STRING_SIZE];

  //erase dummy string content
  memset(string_dummy, '\0', sizeof(string_dummy));
  //copy string content from input string (ignoring blank spaces before and afer string)
  sscanf(string_in,"%s",string_dummy);
  //erase string content
  memset(string_in, '\0', sizeof(&string_in));
  //copy dummy information withou blank spaces in original string
  strcpy(string_in,string_dummy);

  return;
}

void add_object_tolist(char string_name[STRING_SIZE],char string_value[STRING_SIZE], int *number_readobject
                       , char **varname_list,char **value_list) 
{
  //problem: the sscanf line reads between double quotes, including blank space
  //therefore another sscanf is performed in function remove_blankspace_around_string
  //to remove blank spaces in front and after object name (and value resp.)
  remove_blankspaces_around_string(string_name);
  remove_blankspaces_around_string(string_value);
  
  //allocate memory for a new object
  varname_list[*number_readobject] = malloc(STRING_SIZE*sizeof(char *));
  value_list[*number_readobject] = malloc(STRING_SIZE*sizeof(char *));
  //	varname_list[*number_readobject] = (char**)malloc(STRING_SIZE*sizeof(char*));
  //	value_list[*number_readobject] = (char**)malloc(STRING_SIZE*sizeof(char*));

  //copy temp strings into object list
  strcpy(varname_list[*number_readobject],string_name);
  strcpy(value_list[*number_readobject],string_value);

  //count number of read objects
  *number_readobject=*number_readobject+1;
  //log_std("func after number_readobject : %i \n",number_readobject);

  return;
}





