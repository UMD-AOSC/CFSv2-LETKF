# GSW Package
For the conversion to potential temperature from in situ temperature, we use the functions in the TEOS Gibbs-Sea Water (GSW) package:

McDougall, T.J. and P.M. Barker, 2011: Getting started with TEOS-10 and the Gibbs Seawater (GSW)
      Oceanographic Toolbox, 28pp., SCOR/IAPSO WG127, ISBN 978-0-646-55621-5.
            http://www.teos-10.org/software.htm
	    
# Observation Operator

  Penny -
  
  This directory is the initiation of a set of observation operators for each new
  data 'class'. I am assuming that it will be most convenient and efficient to 
  process all similar data type (e.g. profiles, including temperature and salinity
  from Argo, XBTs, CTDs, etc.) in its own observation operator program.

  Ultimately, this should allow distributed work effort on each new data type,
  possibily even allowing the satellite data providers to compute their own
  observation operators for rapid inclusion into the operational ocean data
  assimilation system.
