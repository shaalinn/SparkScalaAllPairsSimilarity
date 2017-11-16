/************************************************************************/
/*! \file da_string.c

\brief Functions for manipulating strings.

Various functions for manipulating strings. Some of these functions 
provide new functionality, whereas others are drop-in replacements
of standard functions (but with enhanced functionality).

Ported from George Karypis' GKlib library by David C. Anastasiu,
with permission, in Aug 2013, then modified.

\date Started 11/1/99
\author George
\version $Id: string.c 14330 2013-05-18 12:15:15Z karypis $
 */
/************************************************************************/

#include "includes.h"



/************************************************************************/
/*! \brief Replaces certain characters in a string.

This function takes a string and replaces all the characters in the
\c fromlist with the corresponding characters from the \c tolist. 
That is, each occurence of <tt>fromlist[i]</tt> is replaced by 
<tt>tolist[i]</tt>. 
If the \c tolist is shorter than \c fromlist, then the corresponding 
characters are deleted. The modifications on \c str are done in place. 
It tries to provide a functionality similar to Perl's \b tr// function.

\param str is the string whose characters will be replaced.
\param fromlist is the set of characters to be replaced.
\param tolist is the set of replacement characters .
\returns A pointer to \c str itself.
 */
/************************************************************************/
char* da_strchr_replace(char* const str, const char* const fromlist, const char* const tolist)
{
	ssize_t i, j, k, len, fromlen, tolen;

	len     = strlen(str);
	fromlen = strlen(fromlist);
	tolen   = strlen(tolist);

	for (i=j=0; i<len; i++) {
		for (k=0; k<fromlen; k++) {
			if (str[i] == fromlist[k]) {
				if (k < tolen)
					str[j++] = tolist[k];
				break;
			}
		}
		if (k == fromlen)
			str[j++] = str[i];
	}
	str[j] = '\0';

	return str;
}




/************************************************************************/
/*! \brief Prunes characters from the end of the string.

This function removes any trailing characters that are included in the
\c rmlist. The trimming stops at the last character (i.e., first character 
from the end) that is not in \c rmlist.  
This function can be used to removed trailing spaces, newlines, etc.
This is a distructive operation as it modifies the string.

\param str is the string that will be trimmed.
\param rmlist contains the set of characters that will be removed.
\returns A pointer to \c str itself.
\sa da_strhprune()
 */
/*************************************************************************/
char* da_strtprune(char* const str, const char* const rmlist)
{
	ssize_t i, j, len;

	len = strlen(rmlist);

	for (i=strlen(str)-1; i>=0; i--) {
		for (j=0; j<len; j++) {
			if (str[i] == rmlist[j])
				break;
		}
		if (j == len)
			break;
	}

	str[i+1] = '\0';

	return str;
}


/************************************************************************/
/*! \brief Prunes characters from the beginning of the string.

This function removes any starting characters that are included in the
\c rmlist. The trimming stops at the first character that is not in 
\c rmlist.
This function can be used to removed leading spaces, tabs, etc.
This is a distructive operation as it modifies the string.

\param str is the string that will be trimmed.
\param rmlist contains the set of characters that will be removed.
\returns A pointer to \c str itself.
\sa da_strtprune()
 */
/*************************************************************************/
char* da_strhprune(char* const str, const char* const rmlist)
{
	ssize_t i, j, len;

	len = strlen(rmlist);

	for (i=0; str[i]; i++) {
		for (j=0; j<len; j++) {
			if (str[i] == rmlist[j])
				break;
		}
		if (j == len)
			break;
	}

	if (i>0) { /* If something needs to be removed */
		for (j=0; str[i]; i++, j++)
			str[j] = str[i];
		str[j] = '\0';
	}

	return str;
}


/************************************************************************/
/*! \brief Converts a string to upper case.

This function converts a string to upper case. This operation modifies the 
string itself.

\param str is the string whose case will be changed.
\returns A pointer to \c str itself.
\sa da_strtolower()
 */
/*************************************************************************/
char* da_strtoupper(char* const str)
{
	int i;

	for (i=0; str[i]!='\0'; str[i]=toupper(str[i]), i++);
	return str;
}


/************************************************************************/
/*! \brief Converts a string to lower case.

This function converts a string to lower case. This operation modifies the 
string itself.

\param str is the string whose case will be changed.
\returns A pointer to \c str itself.
\sa da_strtoupper()
 */
/*************************************************************************/
char* da_strtolower(char* const str)
{
	int i;

	for (i=0; str[i]!='\0'; str[i]=tolower(str[i]), i++);
	return str;
}


/************************************************************************/
/*! \brief Duplicates a string

This function is a replacement for C's standard <em>strdup()</em> function.
The key differences between the two are that da_strdup():
  - uses the dynamic memory allocation routines of \e GKlib. 
  - it correctly handles NULL input strings.

The string that is returned must be freed by da_free().

\param orgstr is the string that will be duplicated.
\returns A pointer to the newly created string.
\sa da_free()
 */
/*************************************************************************/
char* da_strdup(const char* const orgstr)
{
	size_t len;
	char *str=NULL;

	if (orgstr != NULL) {
		len = strlen(orgstr)+1;
		str = (char*) da_malloc(len*sizeof(char), "da_strdup: str");
		strcpy(str, orgstr);
	}

	return str;
}


/************************************************************************/
/*! \brief Case insensitive string comparison.

This function compares two strings for equality by ignoring the case of the
strings. 

\warning This function is \b not equivalent to a case-insensitive 
         <em>strcmp()</em> function, as it does not return ordering 
         information.

\param s1 is the first string to be compared.
\param s2 is the second string to be compared.
\retval 1 if the strings are identical,
\retval 0 otherwise.
 */
/*************************************************************************/
int da_strcasecmp(const char* const s1, const char* const s2)
{
	int i=0;

	if (strlen(s1) != strlen(s2))
		return 0;

	while (s1[i] != '\0') {
		if (tolower(s1[i]) != tolower(s2[i]))
			return 0;
		i++;
	}

	return 1;
}


/************************************************************************/
/*! \brief Compare two strings in revere order

This function is similar to strcmp but it performs the comparison as
if the two strings were reversed.

\param s1 is the first string to be compared.
\param s2 is the second string to be compared.
\retval -1, 0, 1, if the s1 < s2, s1 == s2, or s1 > s2.
 */
/*************************************************************************/
int da_strrcmp(const char* const s1, const char* const s2)
{
	int i1 = strlen(s1)-1;
	int i2 = strlen(s2)-1;

	while ((i1 >= 0) && (i2 >= 0)) {
		if (s1[i1] != s2[i2])
			return (s1[i1] - s2[i2]);
		i1--;
		i2--;
	}

	/* i1 == -1 and/or i2 == -1 */

	if (i1 < i2)
		return -1;
	if (i1 > i2)
		return 1;
	return 0;
}



/************************************************************************/
/*! \brief Converts a time_t time into a string 

This function takes a time_t-specified time and returns a string-formated
representation of the corresponding time. The format of the string is
<em>mm/dd/yyyy hh:mm:ss</em>, in which the hours are in military time.

\param time is the time to be converted.
\return It returns a pointer to a statically allocated string that is 
        over-written in successive calls of this function. If the 
        conversion failed, it returns NULL.

 */
/*************************************************************************/
char* da_time2str(time_t time)
{
	static char datestr[128];
	struct tm *tm;

	tm = localtime(&time);

	if (strftime(datestr, 128, "%m/%d/%Y %H:%M:%S", tm) == 0)
		return NULL;
	else
		return datestr;
}


//
//#if !defined(WIN32) && !defined(__MINGW32__)
///************************************************************************/
///*! \brief Converts a date/time string into its equivalent time_t value
//
//This function takes date and/or time specification and converts it in
//the equivalent time_t representation. The conversion is done using the
//strptime() function. The format that da_str2time() understands is
//<em>mm/dd/yyyy hh:mm:ss</em>, in which the hours are in military time.
//
//\param str is the date/time string to be converted.
//\return If the conversion was successful it returns the time, otherwise
//        it returns -1.
//*/
///*************************************************************************/
//time_t da_str2time(char *str)
//{
//  struct tm time;
//  time_t rtime;
//
//  memset(&time, '\0', sizeof(time));
//
//  if (strptime(str, "%m/%d/%Y %H:%M:%S", &time) == NULL)
//    return -1;
//
//  rtime = mktime(&time);
//  return (rtime < 0 ? 0 : rtime);
//}
//#endif



/*************************************************************************
 * This function returns the key of a particular StringMap ID
 **************************************************************************/
char* da_getStringKey(const da_StringMap_t* const strmap, const char id)
{
	int i;

	for (i=0; strmap[i].name; ++i) {
		if (strmap[i].id == id)
			return strmap[i].name;
	}

	return NULL;
}


/*************************************************************************
 * This function returns the ID of a particular string based on the
 * supplied StringMap array
 **************************************************************************/
int da_getStringID(const da_StringMap_t* const strmap, const char* const key)
{
	int i;

	for (i=0; strmap[i].name; ++i) {
		if (da_strcasecmp(key, strmap[i].name))
			return strmap[i].id;
	}

	return -1;
}
