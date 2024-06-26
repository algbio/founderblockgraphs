/** @file founderblockgraph_cmdline.h
 *  @brief The header file for the command line option parser
 *  generated by GNU Gengetopt version 2.23
 *  http://www.gnu.org/software/gengetopt.
 *  DO NOT modify this file, since it can be overwritten
 *  @author GNU Gengetopt */

#ifndef FOUNDERBLOCKGRAPH_CMDLINE_H
#define FOUNDERBLOCKGRAPH_CMDLINE_H

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h> /* for FILE */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef CMDLINE_PARSER_PACKAGE
/** @brief the program name (used for printing errors) */
#define CMDLINE_PARSER_PACKAGE "founderblockgraphs"
#endif

#ifndef CMDLINE_PARSER_PACKAGE_NAME
/** @brief the complete program name (used for help and version) */
#define CMDLINE_PARSER_PACKAGE_NAME "founderblockgraphs"
#endif

#ifndef CMDLINE_PARSER_VERSION
/** @brief the program version */
#define CMDLINE_PARSER_VERSION "0.5"
#endif

/** @brief Where the command line options are stored */
struct gengetopt_args_info
{
  const char *help_help; /**< @brief Print help and exit help description.  */
  const char *full_help_help; /**< @brief Print help, including hidden options, and exit help description.  */
  const char *version_help; /**< @brief Print version and exit help description.  */
  char * input_arg;	/**< @brief MSA input path.  */
  char * input_orig;	/**< @brief MSA input path original value given at command line.  */
  const char *input_help; /**< @brief MSA input path help description.  */
  char * output_arg;	/**< @brief Index/EFG output path.  */
  char * output_orig;	/**< @brief Index/EFG output path original value given at command line.  */
  const char *output_help; /**< @brief Index/EFG output path help description.  */
  long gap_limit_arg;	/**< @brief Gap limit (suppressed by --elastic) (default='1').  */
  char * gap_limit_orig;	/**< @brief Gap limit (suppressed by --elastic) original value given at command line.  */
  const char *gap_limit_help; /**< @brief Gap limit (suppressed by --elastic) help description.  */
  char * graphviz_output_arg;	/**< @brief Graphviz output path.  */
  char * graphviz_output_orig;	/**< @brief Graphviz output path original value given at command line.  */
  const char *graphviz_output_help; /**< @brief Graphviz output path help description.  */
  char * memory_chart_output_arg;	/**< @brief Memory chart output path.  */
  char * memory_chart_output_orig;	/**< @brief Memory chart output path original value given at command line.  */
  const char *memory_chart_output_help; /**< @brief Memory chart output path help description.  */
  int elastic_flag;	/**< @brief Min-max-length semi-repeat-free segmentation (default=off).  */
  const char *elastic_help; /**< @brief Min-max-length semi-repeat-free segmentation help description.  */
  int gfa_flag;	/**< @brief Saves output in xGFA format (default=off).  */
  const char *gfa_help; /**< @brief Saves output in xGFA format help description.  */
  int output_paths_flag;	/**< @brief Print the original sequences as paths of the xGFA graph (requires --gfa) (default=off).  */
  const char *output_paths_help; /**< @brief Print the original sequences as paths of the xGFA graph (requires --gfa) help description.  */
  char * ignore_chars_arg;	/**< @brief Ignore these characters for the indexability property/pattern matching.  */
  char * ignore_chars_orig;	/**< @brief Ignore these characters for the indexability property/pattern matching original value given at command line.  */
  const char *ignore_chars_help; /**< @brief Ignore these characters for the indexability property/pattern matching help description.  */
  long threads_arg;	/**< @brief Max # threads (default='-1').  */
  char * threads_orig;	/**< @brief Max # threads original value given at command line.  */
  const char *threads_help; /**< @brief Max # threads help description.  */
  long heuristic_subset_arg;	/**< @brief Compute optimal segmentation based on the first ROWNUM MSA rows for performance reasons, then fix the resulting graph iteratively (default='-1').  */
  char * heuristic_subset_orig;	/**< @brief Compute optimal segmentation based on the first ROWNUM MSA rows for performance reasons, then fix the resulting graph iteratively original value given at command line.  */
  const char *heuristic_subset_help; /**< @brief Compute optimal segmentation based on the first ROWNUM MSA rows for performance reasons, then fix the resulting graph iteratively help description.  */
  int disable_elastic_tricks_flag;	/**< @brief Disable the tricks considering the start and end of sequences as unique (default=off).  */
  const char *disable_elastic_tricks_help; /**< @brief Disable the tricks considering the start and end of sequences as unique help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int full_help_given ;	/**< @brief Whether full-help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int input_given ;	/**< @brief Whether input was given.  */
  unsigned int output_given ;	/**< @brief Whether output was given.  */
  unsigned int gap_limit_given ;	/**< @brief Whether gap-limit was given.  */
  unsigned int graphviz_output_given ;	/**< @brief Whether graphviz-output was given.  */
  unsigned int memory_chart_output_given ;	/**< @brief Whether memory-chart-output was given.  */
  unsigned int elastic_given ;	/**< @brief Whether elastic was given.  */
  unsigned int gfa_given ;	/**< @brief Whether gfa was given.  */
  unsigned int output_paths_given ;	/**< @brief Whether output-paths was given.  */
  unsigned int ignore_chars_given ;	/**< @brief Whether ignore-chars was given.  */
  unsigned int threads_given ;	/**< @brief Whether threads was given.  */
  unsigned int heuristic_subset_given ;	/**< @brief Whether heuristic-subset was given.  */
  unsigned int disable_elastic_tricks_given ;	/**< @brief Whether disable-elastic-tricks was given.  */

} ;

/** @brief The additional parameters to pass to parser functions */
struct cmdline_parser_params
{
  int override; /**< @brief whether to override possibly already present options (default 0) */
  int initialize; /**< @brief whether to initialize the option structure gengetopt_args_info (default 1) */
  int check_required; /**< @brief whether to check that all required options were provided (default 1) */
  int check_ambiguity; /**< @brief whether to check for options already specified in the option structure gengetopt_args_info (default 0) */
  int print_errors; /**< @brief whether getopt_long should print an error message for a bad option (default 1) */
} ;

/** @brief the purpose string of the program */
extern const char *gengetopt_args_info_purpose;
/** @brief the usage string of the program */
extern const char *gengetopt_args_info_usage;
/** @brief the description string of the program */
extern const char *gengetopt_args_info_description;
/** @brief all the lines making the help output */
extern const char *gengetopt_args_info_help[];
/** @brief all the lines making the full help output (including hidden options) */
extern const char *gengetopt_args_info_full_help[];

/**
 * The command line parser
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser (int argc, char **argv,
  struct gengetopt_args_info *args_info);

/**
 * The command line parser (version with additional parameters - deprecated)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param override whether to override possibly already present options
 * @param initialize whether to initialize the option structure my_args_info
 * @param check_required whether to check that all required options were provided
 * @return 0 if everything went fine, NON 0 if an error took place
 * @deprecated use cmdline_parser_ext() instead
 */
int cmdline_parser2 (int argc, char **argv,
  struct gengetopt_args_info *args_info,
  int override, int initialize, int check_required);

/**
 * The command line parser (version with additional parameters)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param params additional parameters for the parser
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_ext (int argc, char **argv,
  struct gengetopt_args_info *args_info,
  struct cmdline_parser_params *params);

/**
 * Save the contents of the option struct into an already open FILE stream.
 * @param outfile the stream where to dump options
 * @param args_info the option struct to dump
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_dump(FILE *outfile,
  struct gengetopt_args_info *args_info);

/**
 * Save the contents of the option struct into a (text) file.
 * This file can be read by the config file parser (if generated by gengetopt)
 * @param filename the file where to save
 * @param args_info the option struct to save
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_file_save(const char *filename,
  struct gengetopt_args_info *args_info);

/**
 * Print the help
 */
void cmdline_parser_print_help(void);
/**
 * Print the full help (including hidden options)
 */
void cmdline_parser_print_full_help(void);
/**
 * Print the version
 */
void cmdline_parser_print_version(void);

/**
 * Initializes all the fields a cmdline_parser_params structure 
 * to their default values
 * @param params the structure to initialize
 */
void cmdline_parser_params_init(struct cmdline_parser_params *params);

/**
 * Allocates dynamically a cmdline_parser_params structure and initializes
 * all its fields to their default values
 * @return the created and initialized cmdline_parser_params structure
 */
struct cmdline_parser_params *cmdline_parser_params_create(void);

/**
 * Initializes the passed gengetopt_args_info structure's fields
 * (also set default values for options that have a default)
 * @param args_info the structure to initialize
 */
void cmdline_parser_init (struct gengetopt_args_info *args_info);
/**
 * Deallocates the string fields of the gengetopt_args_info structure
 * (but does not deallocate the structure itself)
 * @param args_info the structure to deallocate
 */
void cmdline_parser_free (struct gengetopt_args_info *args_info);

/**
 * Checks that all the required options were specified
 * @param args_info the structure to check
 * @param prog_name the name of the program that will be used to print
 *   possible errors
 * @return
 */
int cmdline_parser_required (struct gengetopt_args_info *args_info,
  const char *prog_name);


#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* FOUNDERBLOCKGRAPH_CMDLINE_H */
