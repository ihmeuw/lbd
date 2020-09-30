################################################################################
##
## BASE DB TOOLS FOR R
##
## Created: Sept 11, 2018
## Purpose: Functions for making simple connection to MySQL databases,
##   executing arbitrary SQL, and inserting/editing rows using data.frames and
##   data.tables.
##
################################################################################


## Imports for all functions in this section
## NOTE: RMariaDB will be loaded differently in different RStudio images
try(library(RMariaDB), silent=TRUE)
try(library(
  RMariaDB,
  lib.loc="<<<< FILEPATH REDACTED >>>>"
), silent=TRUE)
library(data.table)


## parse_odbc() --------------------------------------------------------------->
#'
#' @title Parse ODBC
#' @description Parses the user's ODBC file and returns details as a nested list
#'
#' @details: Open Database Connectivity (ODBC) is a standard API for accessing
#'   database management systems. This function reads a standard ODBC
#'   initialization file, conventionally named ".odbc.ini", and translates the
#'   results into a named list of database connection definitions.
#'
#' @param odbc_filepath [default='~/.odbc.ini'] Filepath to the ODBC file to be
#'   read.
#'
#' @return Nested list of connection definition arguments
#'
parse_odbc <- function(odbc_filepath='~/.odbc.ini'){
  lines <- scan(odbc_filepath, what='character', sep='\n', quiet=TRUE)
  conn_defs <- list()
  def_name <- ''
  # Iterate through lines, iteratively filling connection definitions
  for (line in lines){
    # Identify whether this is a new connection definition
    name_match <- as.character(grep(pattern='\\[*]', x=line, value=TRUE))
    if (!identical(name_match, character(0))){
      def_name <- substr(x=name_match, start=2, stop=nchar(name_match)-1)
      conn_defs[[def_name]] <- list()
      next
    }
    # Skips any blank leading lines
    if (def_name == ""){
      next
    }
    # Read key, value pairs
    tokens <- strsplit(x=line, split='=')[[1]]
    if (length(tokens)<2){
      next
    }
    k <- tolower(gsub(x=tokens[1],pattern=" ",replace=""))
    v <- gsub(x=tokens[2],pattern=" ",replace="")
    conn_defs[[def_name]][[k]] <- v
  }
  return(conn_defs)
}


## new_db_connection() -------------------------------------------------------->
#'
#' @title New DB connection
#' @description Initialize a new database connection
#'
#' @details: Given ODBC connection definition inputs from \code{\link{parse_odbc}},
#'   create a connection to
#'   a DBMS.
#'
#' @param driver Filepath to an ODBC database driver
#' @param user Character vector giving the DB username
#' @param password Character vector giving the DB password
#' @param server Character vector giving the full path to the DB server
#' @param database Character vector giving the name of the database
#' @param port [default='3306'] The port to connect to on the database
#'
#' @return RMariaDB database connection
#'
#' @seealso \code{\link{parse_odbc}}
#'
new_connection <- function(driver, user, password, server, database,
                           port="<<<< PORT REDACTED >>>>"){
  # Instantiate the RMariaDB connection
  connection <- RMariaDB::dbConnect(drv = RMariaDB::MariaDB(),
                                    host = server,
                                    user = user,
                                    password = password,
                                    dbname = database,
                                    port = port)
  # Return the connection
  return(connection)
}


## connect_from_odbc ---------------------------------------------------------->
#'
#' @title Connect from ODBC
#' @description Connect to a database using an ODBC connection definition
#'
#' @details Given a path to an ODBC initialization file and a "connection
#'   definition" (the name of a connection in that ODBC file), parse the ODBC
#'   file and create a connection to the specified database.

#' @param conn_def The connection definition for a given database
#' @param odbc_filepath [default='~/.odbc.ini'] Filepath to the ODBC file to be
#'   read.
#'
#' @return RMariaDB database connection.
#'
#' @seealso \code{\link{parse_odbc}}, \code{\link{new_db_connection}}
#'
connect_from_odbc <- function(conn_def, odbc_filepath='~/.odbc.ini'){
  # Parse the ODBC file
  credentials <- parse_odbc()[[conn_def]]
  # Create an RMariaDB connection using ODBC file definitions
  connection <- new_connection(driver=credentials[['driver']],
                               user=credentials[['user']],
                               password=credentials[['password']],
                               server=credentials[['server']],
                               database=credentials[['database']],
                               port=credentials[['port']])
  # Return the RMariaDB connection
  return(connection)
}


## execute_sql ---------------------------------------------------------------->
#'
#' @title Execute SQL
#' @description Execute an arbitraray SQL commend in a database
#'
#' @details Given the name of a connection or an ODBC initialization file &
#'   connection definition as well as a SQL command, execute that SQL command
#'   in the specified database.
#'
#' @param statement The SQL statement to execute
#' @param connection [default=NULL] Active RMariaDB connection to the database,
#'   if it exists. Either an active connection OR a valid ODBC filepath &
#'   connection definition must be passed to this function.
#' @param conn_def [default=NULL] Valid connection definition defined within
#'   the ODBC initialization file (set in `odbc_filepath()``). Either a valid
#'   ODBC filepath & connection definition OR an active RMariaDB connection
#'   must be passed to this function.
#' @param odbc_filepath [default='~/.odbc.ini'] Filepath to the ODBC file to be
#'   read. Either a valid ODBC filepath & connection definition OR an active
#'   RMariaDB connection must be passed to this function.
#' @param return_output [default='detect'] Whether or not to return the result
#'   of this function as a data.table. The default, `"detect"`, will return
#'   a data.table for all `"SELECT"` statements and `NULL` otherwise. This
#'   default can be overridden by passing a logical value to this argument
#'
#' @return If `return_output` resolves to `TRUE`, this function will return a
#'   data.table constructed from the results of the database call. If
#'   `return_output` resolves to `FALSE`, this function returns `NULL`.
#'
#' @seealso \code{\link{connect_from_odbc}}
#'
execute_sql <- function(
  statement,
  connection=NULL,
  conn_def=NULL,
  odbc_filepath='~/.odbc.ini',
  return_output='detect'
  ){
  ## Validate inputs
  if( !(return_output %in% c(TRUE, FALSE, 'detect')) ){
    stop("return_output must be one of TRUE, FALSE, or 'detect'.")
  }
  if( !is.null(connection) & (is.null(odbc_filepath)|is.null(conn_def)) ){
    stop(paste0("Either a valid RMariaDB connection or ODBC filepath & ",
                "connection definition must be passed."))
  }
  ## Check whether or not to return the results of the function
  if(return_output=='detect'){
    return_output <- grepl(
      "SELECT",
      toupper(substr(statement,1,15))
    )
  }
  ## Create an RMariaDB connection to the database if one was not passed
  temp_connection <- identical(connection,NULL)
  if (temp_connection==TRUE){
    connection <- connect_from_odbc(
      conn_def=conn_def,
      odbc_filepath=odbc_filepath
    )
  }
  ## Execute SQL statement
  query_result <- RMariaDB::dbSendQuery(connection, statement)
  # Return output, or not
  if (return_output==TRUE){
    query_output <- as.data.table(dbFetch(query_result))
  } else {
    query_output <- NULL
  }
  RMariaDB::dbClearResult(query_result)
  # If the connection was only temporary, close the connection
  if (temp_connection==TRUE){
    RMariaDB::dbDisconnect(connection)
  }
  ## Return output
  return(query_output)
}
