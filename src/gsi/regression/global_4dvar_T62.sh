
set -x

# Set experiment name and analysis date

exp=$jobname

# Set path/file for gsi executable
#gsiexec=$updat

# Set the JCAP resolution which you want.
# All resolutions use LEVS=64
#export JCAP=62
export LEVS=64
export JCAP_B=62

# Set runtime and save directories
tmpdir=$tmpdir/$tmpregdir/${exp}
savdir=$savdir/4dvar_out${JCAP}/sigmap/${exp}

# Specify GSI fixed field and data directories.
fixcrtm=${fixcrtm:-$CRTM_FIX}

# Set variables used in script
#   CLEAN up $tmpdir when finished (YES=remove, NO=leave alone)
#   ncp is cp replacement, currently keep as /bin/cp

UNCOMPRESS=gunzip
CLEAN=NO
ncp=/bin/cp

# Given the requested resolution, set dependent resolution parameters
if [[ "$JCAP" = "382" ]]; then
   export LONA=768
   export LATA=384
   export DELTIM=180
   export resol=1
elif [[ "$JCAP" = "62" ]]; then
   export LONA=192
   export LATA=94
   export DELTIM=1200
   export resol=2
else
   echo "INVALID JCAP = $JCAP"
   exit
fi
export NLAT=$((${LATA}+2))

# Given the analysis date, compute the date from which the
# first guess comes.  Extract cycle and set prefix and suffix
# for guess and observation data files
gdate=`date +%Y%m%d%H -d "${global_4dvar_T62_adate:0:8} ${global_4dvar_T62_adate:8:2} - 6 hours"`
hha=`echo $global_4dvar_T62_adate | cut -c9-10`
hhg=`echo $gdate | cut -c9-10`
prefix_obs=gdas1.t${hha}z
prefix_tbc=gdas1.t${hhg}z
prefix_sfc=gdas${resol}.t${hhg}z
prefix_atm=gdas${resol}.t${hha}z
prefixg=gdas1.t${hhg}z
suffix=tm00.bufr_d

adate0=`echo $global_4dvar_T62_adate | cut -c1-8`
gdate0=`echo $gdate | cut -c1-8`

# Set up $tmpdir
rm -rf $tmpdir
mkdir -p $tmpdir
cd $tmpdir
rm -rf core*

# Make gsi namelist

# CO2 namelist and file decisions
ICO2=${ICO2:-0}
if [ $ICO2 -gt 0 ] ; then
        # Copy co2 files to $tmpdir
        co2dir=${CO2DIR:-$fix_file}
        yyyy=$(echo ${CDATE:-$global_4dvar_T62_adate}|cut -c1-4)
        rm ./global_co2_data.txt
        co2=$co2dir/global_co2.gcmscl_$yyyy.txt
        while [ ! -s $co2 ] ; do
                ((yyyy-=1))
                co2=$co2dir/global_co2.gcmscl_$yyyy.txt
        done
        if [ -s $co2 ] ; then
                $ncp $co2 ./global_co2_data.txt
        fi
        if [ ! -s ./global_co2_data.txt ] ; then
                echo "\./global_co2_data.txt" not created
                exit 1
   fi
fi
#CH4 file decision
ICH4=${ICH4:-0}
if [ $ICH4 -gt 0 ] ; then
#        # Copy ch4 files to $tmpdir
        ch4dir=${CH4DIR:-$fix_file}
        yyyy=$(echo ${CDATE:-$global_4dvar_T62_adate}|cut -c1-4)
        rm ./ch4globaldata.txt
        ch4=$ch4dir/global_ch4_esrlctm_$yyyy.txt
        while [ ! -s $ch4 ] ; do
                ((yyyy-=1))
                ch4=$ch4dir/global_ch4_esrlctm_$yyyy.txt
        done
        if [ -s $ch4 ] ; then
                $ncp $ch4 ./ch4globaldata.txt
        fi
        if [ ! -s ./ch4globaldata.txt ] ; then
                echo "\./ch4globaldata.txt" not created
                exit 1
   fi
fi
IN2O=${IN2O:-0}
if [ $IN2O -gt 0 ] ; then
#        # Copy ch4 files to $tmpdir
        n2odir=${N2ODIR:-$fix_file}
        yyyy=$(echo ${CDATE:-$global_4dvar_T62_adate}|cut -c1-4)
        rm ./n2oglobaldata.txt
        n2o=$n2odir/global_n2o_esrlctm_$yyyy.txt
        while [ ! -s $n2o ] ; do
                ((yyyy-=1))
                n2o=$n2odir/global_n2o_esrlctm_$yyyy.txt
        done
        if [ -s $n2o ] ; then
                $ncp $n2o ./n2oglobaldata.txt
        fi
        if [ ! -s ./n2oglobaldata.txt ] ; then
                echo "\./n2oglobaldata.txt" not created
                exit 1
   fi
fi
ICO=${ICO:-0}
if [ $ICO -gt 0 ] ; then
#        # Copy CO files to $tmpdir
        codir=${CODIR:-$fix_file}
        yyyy=$(echo ${CDATE:-$global_4dvar_T62_adate}|cut -c1-4)
        rm ./coglobaldata.txt
        co=$codir/global_co_esrlctm_$yyyy.txt
        while [ ! -s $co ] ; do
                ((yyyy-=1))
                co=$codir/global_co_esrlctm_$yyyy.txt
        done
        if [ -s $co ] ; then
                $ncp $co ./coglobaldata.txt
        fi
        if [ ! -s ./coglobaldata.txt ] ; then
                echo "\./coglobaldata.txt" not created
                exit 1
   fi
fi

. $scripts/regression_nl_update.sh

GRIDOPTS="$GRIDOPTS_update"
BKGVERR="$BKGVERR_update"
ANBKGERR="$ANBKERR_update"
JCOPTS="$JCOPTS_update"
STRONGOPTS="$STRONGOPTS_update"
OBSQC="$OBSQC_update"
OBSINPUT="$OBSINPUT_update"
SUPERRAD="$SUPERRAD_update"
SINGLEOB="$SINGLEOB_update"

# Set variables for requested minimization (pcgsoi or lanczos)
JCOPTS="ljcpdry=.false.,"
OBSQC="noiqc=.false.,"
SETUPmin="miter=1,niter(1)=50,niter_no_qc(1)=500,"
SETUPlan=""
#export minimization=${minimization:-"pcgsoi"}
#if [ "$minimization" = "lanczos" ]; then
#   SETUPlan="lsqrtb=.true.,lcongrad=.true.,ltlint=.true.,ladtest=.true.,lgrtest=.false.,"
#   HYBENS_GLOBAL=".false."
#fi

# Create namelist for observer run
export nhr_obsbin=${nhr_obsbin:-1}
SETUPobs="l4dvar=.true.,jiterstart=1,lobserver=.true.,iwrtinc=1,nhr_assimilation=6,nhr_obsbin=$nhr_obsbin,"
#SETUP="$SETUPmin $SETUPlan $SETUPobs $SETUP_update"
SETUP="$SETUPmin $SETUPobs $SETUP_update"

if [ "$minimization" = "lanczos" ]; then
   namelist_name=global_lanczos_T62
else
   namelist_name=global_T62
fi

if [ "$debug" = ".false." ]; then
   . $scripts/regression_namelists.sh $namelist_name
else
   . $scripts/regression_namelists_db.sh $namelist_name
fi
rm gsiparm.anl
cat << EOF > gsiparm.anl

$gsi_namelist

EOF


# Set fixed files
#   berror   = forecast model background error statistics
#   specoef  = CRTM spectral coefficients
#   trncoef  = CRTM transmittance coefficients
#   emiscoef = CRTM coefficients for IR sea surface emissivity model
#   aerocoef = CRTM coefficients for aerosol effects
#   cldcoef  = CRTM coefficients for cloud effects
#   satinfo  = text file with information about assimilation of brightness temperatures
#   cloudyinfo  = text file with information about assimilation of cloudy radiance
#   satangl  = angle dependent bias correction file (fixed in time)
#   atmsbeamdat  =  data required for atms spatial averaging
#   pcpinfo  = text file with information about assimilation of prepcipitation rates
#   ozinfo   = text file with information about assimilation of ozone data
#   errtable = text file with obs error for conventional data (optional)
#   convinfo = text file with information about assimilation of conventional data
#   bufrtable= text file ONLY needed for single obs test (oneobstest=.true.)
#   bftab_sst= bufr table for sst ONLY needed for sst retrieval (retrieval=.true.)

anavinfo=$fixgsi/global_anavinfo.l64.txt
berror=$fixgsi/$endianness/global_berror.l${LEVS}y${NLAT}.f77
emiscoef_IRwater=$fixcrtm/Nalli.IRwater.EmisCoeff.bin
emiscoef_IRice=$fixcrtm/NPOESS.IRice.EmisCoeff.bin
emiscoef_IRland=$fixcrtm/NPOESS.IRland.EmisCoeff.bin
emiscoef_IRsnow=$fixcrtm/NPOESS.IRsnow.EmisCoeff.bin
emiscoef_VISice=$fixcrtm/NPOESS.VISice.EmisCoeff.bin
emiscoef_VISland=$fixcrtm/NPOESS.VISland.EmisCoeff.bin
emiscoef_VISsnow=$fixcrtm/NPOESS.VISsnow.EmisCoeff.bin
emiscoef_VISwater=$fixcrtm/NPOESS.VISwater.EmisCoeff.bin
emiscoef_MWwater=$fixcrtm/FASTEM6.MWwater.EmisCoeff.bin
aercoef=$fixcrtm/AerosolCoeff.bin
cldcoef=$fixcrtm/CloudCoeff.bin
satinfo=$fixgsi/global_satinfo.txt
cloudyinfo=$fixgsi/cloudy_radiance_info.txt
scaninfo=$fixgsi/global_scaninfo.txt
satangl=$fixgsi/global_satangbias.txt
atmsbeamdat=$fixgsi/atms_beamwidth.txt
pcpinfo=$fixgsi/global_pcpinfo.txt
ozinfo=$fixgsi/global_ozinfo.txt
convinfo=$fixgsi/global_convinfo_reg_test.txt
vqcdat=$fixgsi/vqctp001.dat
errtable=$fixgsi/prepobs_errtable.global
### add 9 tables
errtable_pw=$fixgsi/prepobs_errtable_pw.global
errtable_ps=$fixgsi/prepobs_errtable_ps.global_nqcf
errtable_t=$fixgsi/prepobs_errtable_t.global_nqcf
errtable_q=$fixgsi/prepobs_errtable_q.global_nqcf
errtable_uv=$fixgsi/prepobs_errtable_uv.global_nqcf
btable_ps=$fixgsi/nqc_b_ps.global_nqcf
btable_t=$fixgsi/nqc_b_t.global_nqcf
btable_q=$fixgsi/nqc_b_q.global_nqcf
btable_uv=$fixgsi/nqc_b_uv.global_nqcf


# Only need this file for single obs test
bufrtable=$fixgsi/prepobs_prep.bufrtable

# Only need this file for sst retrieval
bftab_sst=$fixgsi/bufrtab.012

# Copy executable and fixed files to $tmpdir
if [[ $exp == *"updat"* ]]; then
   $ncp $gsiexec_updat  ./gsi.x
elif [[ $exp == *"contrl"* ]]; then
   $ncp $gsiexec_contrl ./gsi.x
fi

$ncp $anavinfo ./anavinfo
$ncp $berror   ./berror_stats
$ncp $emiscoef_IRwater ./Nalli.IRwater.EmisCoeff.bin
$ncp $emiscoef_IRice ./NPOESS.IRice.EmisCoeff.bin
$ncp $emiscoef_IRsnow ./NPOESS.IRsnow.EmisCoeff.bin
$ncp $emiscoef_IRland ./NPOESS.IRland.EmisCoeff.bin
$ncp $emiscoef_VISice ./NPOESS.VISice.EmisCoeff.bin
$ncp $emiscoef_VISland ./NPOESS.VISland.EmisCoeff.bin
$ncp $emiscoef_VISsnow ./NPOESS.VISsnow.EmisCoeff.bin
$ncp $emiscoef_VISwater ./NPOESS.VISwater.EmisCoeff.bin
$ncp $emiscoef_MWwater ./FASTEM6.MWwater.EmisCoeff.bin
$ncp $aercoef  ./AerosolCoeff.bin
$ncp $cldcoef  ./CloudCoeff.bin
$ncp $satangl  ./satbias_angle
$ncp $atmsbeamdat  ./atms_beamwidth.txt
$ncp $satinfo  ./satinfo
$ncp $cloudyinfo  ./cloudy_radiance_info.txt
$ncp $scaninfo ./scaninfo
$ncp $pcpinfo  ./pcpinfo
$ncp $ozinfo   ./ozinfo
$ncp $convinfo ./convinfo
$ncp $vqcdat ./vqctp001.dat
$ncp $errtable ./errtable
#add 9 tables for new varqc
$ncp $errtable_pw           ./errtable_pw
$ncp $errtable_ps           ./errtable_ps
$ncp $errtable_t           ./errtable_t
$ncp $errtable_q           ./errtable_q
$ncp $errtable_uv           ./errtable_uv
$ncp $btable_ps           ./btable_ps
$ncp $btable_t           ./btable_t
$ncp $btable_q           ./btable_q
$ncp $btable_uv           ./btable_uv



$ncp $bufrtable ./prepobs_prep.bufrtable
$ncp $bftab_sst ./bftab_sstphr

#if using correlated error, link to the covariance files
#if grep -q "Rcov" $anavinfo ;
#then 
#  if ls ${fixgsi}/Rcov* 1> /dev/null 2>&1;
#  then
#    $ncp ${fixgsi}/Rcov* .
#  else
#    echo "Warning: Satellite error covariance files are missing."
#    echo "Check for the required Rcov files in " $anavinfo
#    exit 1
#  fi
#fi

# Adjust data usage flags in convinfo file.
rm new
cp convinfo old
mv convinfo convinfo_original
sed 's/sst      180    0   -1     3.0/sst      180    0    1     3.0/' < old > new
mv new old
sed 's/uv       243   56    1     3.0/uv       243   56   -1     3.0/' < old > new
mv new old
sed 's/uv       253   56    1     3.0/uv       253   56   -1     3.0/' < old > new
mv new convinfo

# Copy CRTM coefficient files based on entries in satinfo file
for file in `awk '{if($1!~"!"){print $1}}' ./satinfo | sort | uniq` ;do
   $ncp $fixcrtm/${file}.SpcCoeff.bin ./
   $ncp $fixcrtm/${file}.TauCoeff.bin ./
done

# Copy observational data to $tmpdir
$ncp $global_4dvar_T62_obs/${prefix_obs}.prepbufr                ./prepbufr
$ncp $global_4dvar_T62_obs/${prefix_obs}.satwnd.${suffix}        ./satwndbufr
$ncp $global_4dvar_T62_obs/${prefix_obs}.gpsro.${suffix}         ./gpsrobufr
$ncp $global_4dvar_T62_obs/${prefix_obs}.spssmi.${suffix}        ./ssmirrbufr
$ncp $global_4dvar_T62_obs/${prefix_obs}.sptrmm.${suffix}        ./tmirrbufr
$ncp $global_4dvar_T62_obs/${prefix_obs}.osbuv8.${suffix}        ./sbuvbufr
$ncp $global_4dvar_T62_obs/${prefix_obs}.goesfv.${suffix}        ./gsnd1bufr
$ncp $global_4dvar_T62_obs/${prefix_obs}.1bamua.${suffix}        ./amsuabufr
$ncp $global_4dvar_T62_obs/${prefix_obs}.1bamub.${suffix}        ./amsubbufr
$ncp $global_4dvar_T62_obs/${prefix_obs}.1bhrs2.${suffix}        ./hirs2bufr
$ncp $global_4dvar_T62_obs/${prefix_obs}.1bhrs3.${suffix}        ./hirs3bufr
$ncp $global_4dvar_T62_obs/${prefix_obs}.1bhrs4.${suffix}        ./hirs4bufr
$ncp $global_4dvar_T62_obs/${prefix_obs}.1bmhs.${suffix}         ./mhsbufr
$ncp $global_4dvar_T62_obs/${prefix_obs}.1bmsu.${suffix}         ./msubufr
$ncp $global_4dvar_T62_obs/${prefix_obs}.airsev.${suffix}        ./airsbufr
$ncp $global_4dvar_T62_obs/${prefix_obs}.sevcsr.${suffix}        ./seviribufr
$ncp $global_4dvar_T62_obs/${prefix_obs}.mtiasi.${suffix}        ./iasibufr
$ncp $global_4dvar_T62_obs/${prefix_obs}.ssmit.${suffix}         ./ssmitbufr
$ncp $global_4dvar_T62_obs/${prefix_obs}.amsre.${suffix}         ./amsrebufr
$ncp $global_4dvar_T62_obs/${prefix_obs}.ssmis.${suffix}         ./ssmisbufr
$ncp $global_4dvar_T62_obs/${prefix_obs}.gome.${suffix}          ./gomebufr
$ncp $global_4dvar_T62_obs/${prefix_obs}.omi.${suffix}           ./omibufr
$ncp $global_4dvar_T62_obs/${prefix_obs}.mlsbufr.${suffix}       ./mlsbufr
$ncp $global_4dvar_T62_obs/${prefix_obs}.eshrs3.${suffix}        ./hirs3bufrears
$ncp $global_4dvar_T62_obs/${prefix_obs}.esamua.${suffix}        ./amsuabufrears
$ncp $global_4dvar_T62_obs/${prefix_obs}.esamub.${suffix}        ./amsubbufrears
$ncp $global_4dvar_T62_obs/${prefix_obs}.syndata.tcvitals.tm00   ./tcvitl

# Copy bias correction, atmospheric and surface files
if [ "$minimization" = "lanczos" ]; then
   $ncp $global_4dvar_T62_ges/${prefix_tbc}.abias.orig        ./satbias_in
   $ncp $global_4dvar_T62_ges/${prefix_tbc}.satang.orig       ./satbias_angle
else
   $ncp $global_4dvar_T62_ges/${prefix_tbc}.abias             ./satbias_in
   $ncp $global_4dvar_T62_ges/${prefix_tbc}.abias_pc          ./satbias_pc
   $ncp $global_4dvar_T62_ges/${prefix_tbc}.satang            ./satbias_angle
   $ncp $global_4dvar_T62_ges/${prefix_tbc}.radstat           ./radstat.gdas

   listdiag=`tar xvf radstat.gdas | cut -d' ' -f2 | grep _ges`
   for type in $listdiag; do
      diag_file=`echo $type | cut -d',' -f1`
      fname=`echo $diag_file | cut -d'.' -f1`
      date=`echo $diag_file | cut -d'.' -f2`
      $UNCOMPRESS $diag_file
      fnameanl=$(echo $fname|sed 's/_ges//g')
      mv $fname.$date $fnameanl
   done
fi

if [[ "$endianness" = "Big_Endian" ]]; then
   ##$ncp $global_4dvar_T62_ges/${prefix_sfc}.bf03               ./sfcf03
   $ncp $global_4dvar_T62_ges/${prefix_sfc}.bf06                 ./sfcf06
   ##$ncp $global_4dvar_T62_ges/${prefix_sfc}.bf09               ./sfcf09
elif [[ "$endianness" = "Little_Endian" ]]; then
   ##$ncp $global_4dvar_T62_ges/${prefix_sfc}.bf03.le            ./sfcf03
   $ncp $global_4dvar_T62_ges/${prefix_sfc}.bf06.le              ./sfcf06
   ##$ncp $global_4dvar_T62_ges/${prefix_sfc}.bf09.le            ./sfcf09
fi

if [[ "$endianness" = "Big_Endian" ]]; then
   ##$ncp $global_4dvar_T62_obs/${prefix_atm}.sgm3prep           ./sigf03
   $ncp $global_4dvar_T62_obs/${prefix_atm}.sgesprep             ./sigf06
   ##$ncp $global_4dvar_T62_obs/${prefix_atm}.sgp3prep           ./sigf09
elif [[ "$endianness" = "Little_Endian" ]]; then
   ##$ncp $global_4dvar_T62_obs/${prefix_atm}.sgm3prep.le        ./sigf03
   $ncp $global_4dvar_T62_obs/${prefix_atm}.sgesprep.le          ./sigf06
   ##$ncp $global_4dvar_T62_obs/${prefix_atm}.sgp3prep.le        ./sigf09
fi

# Run gsi observer
cd $tmpdir
echo "run gsi now"
eval "$APRUN $tmpdir/gsi.x > stdout.obsvr 2>&1"

# Run gsi identity model 4dvar under Parallel Operating Environment (poe) on NCEP IBM
rm -f siganl sfcanl.gsi satbias_out fort.2*
rm -rf dir.0*

# Create namelist for identity model 4dvar run
SETUP4dv="l4dvar=.true.,jiterstart=1,nhr_assimilation=6,nhr_obsbin=$nhr_obsbin,idmodel=.true.,iwrtinc=1,lanczosave=.true.,"
#SETUP="$SETUPmin $SETUPlan $SETUP4dv $SETUP_update"
SETUP="$SETUPmin $SETUP4dv $SETUP_update"

if [ "$minimization" = "lanczos" ]; then
   namelist_name=global_lanczos_T62
else
   namelist_name=global_T62
fi

if [ "$debug" = ".false." ]; then
   . $scripts/regression_namelists.sh $namelist_name
else
   . $scripts/regression_namelists_db.sh $namelist_name
fi
rm gsiparm.anl
cat << EOF > gsiparm.anl
$gsi_namelist
EOF

cd $tmpdir
echo "run gsi now"
eval "$APRUN $tmpdir/gsi.x > stdout 2>&1"
rc=$?
exit $rc
