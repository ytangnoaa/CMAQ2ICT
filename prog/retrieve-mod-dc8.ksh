#!/bin/ksh

set -x 
set -A fdate 20190722 20190724 20190725 20190729 20190730 20190802 20190803 \
      20190806 20190807 20190808 20190812 20190813 20190815 20190816 20190819 \
      20190821 20190823 20190826 20190829 20190830 20190831 20190903 20190905

set -A fcase 1 1

sflight=11
eflight=11
#sflight=22
#eflight=23

runextract=1
runtest=1

HOLDIN1=$HOME/ptmp2/com/aqm/para/aqm.
HOLDIN2=$HOME/ptmp2/com/aqm/wrf/aqm.

HPSS1=/5year/NCEPDEV/emc-naqfc/Youhua.Tang/cmaq531-test/firevoc.
HPSS2=/5year/NCEPDEV/emc-naqfc/Youhua.Tang/cmaq531-test/wrfcmaq.

export IOAPI_ISPH=20
export IOAPI_CHECK_HEADERS=F

BASE=$(pwd)

function checkfile
{
 nowfile=dc8-1m-${sflight}.dat

 firex-dc8-5x-cb6.x $sflight
aline=`wc -c $nowfile`
if [ $? -gt 0 ]; then
  echo "PARA $nowfile does not exist"
  exit 2
fi
echo $aline
for var in $aline ; do
  asize=$var
  break;
done
if [ $asize -lt 400000 ]; then
  echo "PARA $nowfile wrong size $asize "
  exit 2
fi
}


while [ $sflight -le $eflight ]; do

export GRIDDESC=/gpfs/hps3/emc/naqfc/noscrub/Youhua.Tang/nwdev/NAQFC-WCOSS/parm/aqm_griddesc05
export TOPO=aqm.t12z.grdcro2d.ncf
export TOPODOT=aqm.t12z.grddot2d.ncf


 nowdate=${fdate[sflight-1]}
 nowdatem1=$(ndate -24 ${nowdate}00 |cut -c1-8)
 nowyear=$(echo $nowdate|cut -c1-4)
 echo "nowdate = $nowdate $sflight"

if [ ${fcase[0]} -eq 1 ]; then

if [ ! -s $HOLDIN1$nowdate/aqm.t12z.conc.ncf ]; then
mkdir -p $HOLDIN1$nowdate

hsi <<EOF
cd $HPSS1$nowdate
lcd $HOLDIN1$nowdate
get aqm.t12z.conc.ncf
get aqm.t12z.rj_2.ncf
get aqm.t12z.rj_3.ncf
get aqm.t12z.pmdiag.ncf
get aqm.t12z.met*3d.ncf
bye
EOF
fi
if [ ! -s $HOLDIN1$nowdate/aqm.t12z.metcro3d.ncf ]; then
  cd $HOLDIN1$nowdate
  hpsstar get /NAGAPE/arl/5year/Patrick.C.Campbell/NACC_FV3GFS16_runs/base_nacc_cmaq531_nofire/GFSv16.NACCv100.${nowdate}12.output_met.tar
   if [ -s $HOLDIN1$nowdate/aqm.t12z.metcro3d.ncf ]; then
    hsi "cd $HPSS1$nowdate; put aqm.t12z.met*ncf ; put aqm*soi*ncf; put aqm*lufraccro*ncf" &
   else
    echo " can not metcro3d"
    exit 1 
   fi
fi

export WINDDOT=$HOLDIN1$nowdate/aqm.t12z.metdot3d.ncf
export WINDDOT2=$WINDDOT
export MET3DCRO=$HOLDIN1$nowdate/aqm.t12z.metcro3d.ncf
export MET3DCRO2=$MET3DCRO
export CHEM3D=$HOLDIN1$nowdate/aqm.t12z.conc.ncf
export CHEM3D2=$CHEM3D
export PMDIAG=$HOLDIN1$nowdate/aqm.t12z.pmdiag.ncf
export PMDIAG2=$PMDIAG
export JVFILE=$HOLDIN1$nowdate/aqm.t12z.rj_2.ncf
export JVFILE2=$JVFILE
export AOP=$HOLDIN1$nowdate/aqm.t12z.rj_3.ncf
export AOP2=$AOP

cd $BASE
typeset -Z2 sflight
checkfile
mv $nowfile ../data/dc8-1m-firevoc-${sflight}.dat

fi

# WRF-CMAQ
if [ ${fcase[1]} -eq 1 ]; then

if [ ! -s $HOLDIN2$nowdate/aqm.t12z.conc.ncf ]; then
mkdir -p $HOLDIN2$nowdate 
hsi <<EOF
prompt
cd $HPSS2$nowdate
lcd $HOLDIN2$nowdate
get aqm.t12z.conc.ncf
get aqm.t12z.rj_2.ncf
get aqm.t12z.rj_3.ncf
get aqm.t12z.pmdiag.ncf
get aqm.t12z.met*3d.ncf
bye
EOF
fi
if [ ! -s $HOLDIN2$nowdate/aqm.t12z.conc.ncf ] || [ ! -s $HOLDIN2$nowdate/aqm.t12z.metcro3d.ncf ] ; then
 echo " miss files "
 exit 1
fi 
export WINDDOT=$HOLDIN2$nowdate/aqm.t12z.metdot3d.ncf
export WINDDOT2=$WINDDOT
export MET3DCRO=$HOLDIN2$nowdate/aqm.t12z.metcro3d.ncf
export MET3DCRO2=$MET3DCRO
export CHEM3D=$HOLDIN2$nowdate/aqm.t12z.conc.ncf
export CHEM3D2=$CHEM3D
export PMDIAG=$HOLDIN2$nowdate/aqm.t12z.pmdiag.ncf
export PMDIAG2=$PMDIAG
export JVFILE=$HOLDIN2$nowdate/aqm.t12z.rj_2.ncf
export JVFILE2=$JVFILE
export AOP=$HOLDIN2$nowdate/aqm.t12z.rj_3.ncf
export AOP2=$AOP

cd $BASE
typeset -Z2 sflight
checkfile
mv $nowfile ../data/dc8-1m-wrfcmaq-${sflight}.dat
fi

 let sflight=sflight+1
done
