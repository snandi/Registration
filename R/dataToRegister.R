# setClass for "dataToRegister"

dataToRegister <- setClass( 
  # Set the name for the class
  Class = "dataToRegister", 

  # Define the slots
  slots = c(
    curves     = "matrix",   ## a matrix or dataframe with each curve in a separate column
    abscissa   = "numeric",  ## abscissa, or x-axis, or pixels
    curveNames = "character" ## column names
  ),
  
  # Make a function that can test to see if the data is consistent.
  # This is not called if you have an initialize function defined!
  validity = function( object )
  {
    if( class( object@curves ) != "matrix" | class( object@curves ) != "data.frame" ){
      return("Curves to register must be of matrix or data.frame class")
    }
    if( class( abscissa ) != "numeric" ){
      return("abscissa must be a numeric vector")
    }
    if( class( curveNames ) != "character" ){
      return("abscissa must be a character vector")
    }
    return(TRUE)
  }
  
)
