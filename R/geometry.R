# calcChamberGeometry calcualtes inner chamber volume and enclosed respiring surface area given the geometry and dimensions of the chamber.
# To do: account for curbed respiring surface, such as tree stems

calcChamberGeometry <- function (
  geometry   = 'cylinder', # Geometry of the chamber such as cylinder or cuboid
  dimensions = NA,         # Vector of chamber dimensions in meters
                           # For a cylinder, first element = radius, second element = height 
                           # For a cuboid,   first element = width,  second element = depth, third element = height
  taper      = 1.0         # A form fctor describing a linear taper of the chamber. 1.0 is equal to no taper. 
){
  if (geometry == 'cylinder') {
    RespArea      <- pi * dimensions [1]**2.0
    ChamberVolume <- RespArea * dimensions [2] * taper 
  } else if (geometry == 'cuboid') {
    RespArea      <- dimensions [1] * dimensions [2]
    ChamberVolume <- RespArea * dimensions [3] * taper 
  }
  c (ChamberVolume, RespArea)
}

attr (calcChamberGeometry,"ex") <- function(){
  innerRadius <- 0.1016
  innerHeight <- 0.0762
  calcChamberGeometry (geometry = 'cylinder', dimensions = c (innerRadius, innerHeight))
}
