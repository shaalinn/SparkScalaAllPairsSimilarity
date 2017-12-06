/*!
\file da_getopt.h
\brief This file contains GNU's externs/structs/prototypes

Ported from GKlib, by George Karypis

\date   Started 3/27/2007
\author George
\version\verbatim $Id: da_getopt.h 10711 2011-08-31 22:23:04Z karypis $ \endverbatim
 */

#ifndef _COSGRAPH_GETOPT_H_
#define _COSGRAPH_GETOPT_H_


/* Externals from getopt.c */
extern char *da_optarg;
extern int da_optind;
extern int da_opterr;
extern int da_optopt;


/*! \brief The structure that stores the information about the command-line options 

This structure describes a single long option name for the sake of 
da_getopt_long(). The argument <tt>long_options</tt> must be an array
of these structures, one for each long option. Terminate the array with 
an element containing all zeros.
 */
struct da_option
{
  char *name;			/*!< This field is the name of the option. */
  int has_arg;			/*!< This field says whether the option takes an argument.
				   It is an integer, and there are three legitimate values: 
				   no_argument, required_argument and optional_argument. 
				 */
  int *flag;			/*!< See the discussion on ::da_option#val */
  int val;			/*!< These fields control how to report or act on the option
				   when it occurs. 

				   If flag is a null pointer, then the val is a value which 
				   identifies this option. Often these values are chosen 
				   to uniquely identify particular long options.

				   If flag is not a null pointer, it should be the address 
				   of an int variable which is the flag for this option. 
				   The value in val is the value to store in the flag to 
				   indicate that the option was seen. */
};

/* Names for the values of the `has_arg' field of `struct da_option'.  */
#define no_argument		0
#define required_argument	1
#define optional_argument	2


/* Function prototypes */
extern int da_getopt (int __argc, char **__argv, char *__shortopts);
extern int da_getopt_long (int __argc, char **__argv, char *__shortopts,
			   struct da_option *__longopts, int *__longind);
extern int da_getopt_long_only (int __argc, char **__argv,
				char *__shortopts,
				struct da_option *__longopts, int *__longind);

#endif
