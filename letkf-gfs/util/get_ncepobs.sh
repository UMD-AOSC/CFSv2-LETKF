#!/bin/bash
set -e
#===============================================================================
#
#  Download NCEP conventional observation data from UCAR/DSS
#   -- adapted from Takemasa Miyoshi's LETKF google code,
#      March 2013, Guo-Yuan Lien
#
#===============================================================================
# if [ -f configure.sh ]; then
#   . configure.sh
# else
#   echo "[Error] $0: 'configure.sh' does not exist." 1>&2
#   exit 1
# fi

if [ "$#" -lt 3 ]; then
  cat 1>&2 << EOF

[get_ncepobs.sh] Download NCEP conventional observation data.
                 *use settings in 'configure.sh'

Usage: $0 EMAIL PASSWD STIME [ETIME] [IF_DECODE]

  EMAIL      UCAR/DSS account (register at http://rda.ucar.edu )
  PASSWD     UCAR/DSS account password
  STIME      Start time (format: YYYYMMDD)
  ETIME      End   time (format: YYYYMMDD)
             (default: same as STIME)
  IF_DECODE  Convert PREPBUFR to LETKF obs format or not
             0: No,  store only PREPBUFR format
             1: Yes, decode to LETKF obs format
             (default: Yes)

EOF
  exit 1
fi

EMAIL="$1"
PASSWD="$2"
STIME=$3
if [ "${#STIME}" -ne "8" ]; then
    echo "STIME needs to be in YYYYMMDD format"
    exit 1
fi
ETIME=${4:-$STIME}
if [ "${#ETIME}" -ne "8" ]; then
    echo "ETIME needs to be in YYYYMMDD format"
    exit 1
fi
IF_DECODE=${5:-1}

BUFRBIN=$CFS_LETKF_ROOT/letkf-gfs/util
DATAURL="http://rda.ucar.edu/data/ds337.0/tarfiles"
LOGINURL="https://rda.ucar.edu/cgi-bin/login"
WGET="/usr/bin/wget --no-check-certificate"
AUTH="auth.rda_ucar_edu"
OPTLOGIN="-O /dev/null --save-cookies $AUTH --post-data=\"email=${EMAIL}&passwd=${PASSWD}&action=login\""
OPT="-N --load-cookies $AUTH"

tmpsubdir="gfs-letkf_${USER}_${SYSNAME}_get_ncepobs"
tmprun="$TMP_DIR_LOCAL/${tmpsubdir}"

#===============================================================================

mkdir -p $tmprun
rm -fr $tmprun/*
mkdir -p $tmprun/download
cd $tmprun
if [ "$IF_DECODE" = '1' ]; then
  cp $BUFRBIN/dec_prepbufr .
fi

cd $tmprun/download
$WGET $OPTLOGIN $LOGINURL

time=$STIME
while [ "$time" -le "$ETIME" ]; do

  yyyy=${time:0:4}
  mm=${time:4:2}
  dd=${time:6:2}
  if [ "$yyyy$mm" -ge 200807 ]; then
    DATAF="prepbufr.$yyyy$mm$dd.nr.tar.gz"
  else
    DATAF="prepbufr.$yyyy$mm$dd.wo40.tar.gz"
  fi

  cd $tmprun/download
  $WGET $OPT "${DATAURL}/$yyyy/${DATAF}"
  tar xzf $DATAF
  rm -f $DATAF

  cd $tmprun
  for hh in '00' '06' '12' '18'; do
    timef="$yyyy$mm$dd$hh"
    echo
    echo "[${timef}]"
    echo

    if [ "$yyyy$mm" -ge 200807 ]; then
      mv -f download/$yyyy$mm$dd.nr/prepbufr.gdas.$yyyy$mm$dd.t${hh}z.nr \
            prepbufr.gdas.${timef}.nr
    else
      mv -f download/$yyyy$mm$dd.wo40/prepbufr.gdas.${timef}.wo40 \
            prepbufr.gdas.${timef}.nr
    fi
    wc -c prepbufr.gdas.${timef}.nr | $BUFRBIN/grabbufr prepbufr.gdas.${timef}.nr prepbufr.in

    if [ "$IF_DECODE" = '1' ]; then
      time $BUFRBIN/dec_prepbufr
      touch fort.87
      touch fort.88
      touch fort.89
      touch fort.90
      touch fort.91
      touch fort.92
      touch fort.93
      obsdir=$OBS_ATM/$yyyy/$yyyy$mm/$yyyy$mm$dd/${timef}
      mkdir -p $obsdir
      mv fort.87 $obsdir/t-3.dat
      mv fort.88 $obsdir/t-2.dat
      mv fort.89 $obsdir/t-1.dat
      mv fort.90 $obsdir/t.dat
      mv fort.91 $obsdir/t+1.dat
      mv fort.92 $obsdir/t+2.dat
      mv fort.93 $obsdir/t+3.dat

    else
        mkdir -p $OBSNCEP/obs${timef}
	mv -f prepbufr.gdas.${timef}.nr $OBSNCEP/obs${timef}
    fi
    
#    mv -f prepbufr.gdas.${timef}.nr $OBSNCEP/obs${timef}/gdas1.t${hh}z.prepbufr.nr
  done
  rm $tmprun/prepbufr.gdas.*.nr
  rm $tmprun/download/*.nr -r

  
time=$(date +%Y%m%d -d "$time + 1 day")
done

rm -rf $tmprun

#===============================================================================

exit 0
