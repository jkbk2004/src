
set -x

# Set experiment name and analysis date

exp=$jobname

# Set path/file for gsi executable
#basedir=/scratch1/portfolios/NCEPDEV/da/save/Daryl.Kleist
#gsipath=$basedir/gsi/
#gsiexec=$gsipath/trunk/src/global_gsi.x

# Set the JCAP resolution which you want.
# All resolutions use LEVS=64
export JCAP=126
export LEVS=64
export JCAP_B=126
export JCAP_EN=62

# Set runtime and save directories
tmpdir=$tmpdir/$tmpregdir/${exp}
savdir=$savdir/out${JCAP}/${exp}

# Specify GSI fixed field and data directories.
#fixcrtm=${fixcrtm:-$CRTM_FIX}

#datobs=/scratch1/portfolios/NCEPDEV/da/noscrub/Daryl.Kleist/CASES/$adate/obs
#datges=/scratch1/portfolios/NCEPDEV/da/noscrub/Daryl.Kleist/CASES/$adate/ges
#datens=/scratch1/portfolios/NCEPDEV/da/noscrub/Daryl.Kleist/CASES/$adate/ens

# Set variables used in script
#   CLEAN up $tmpdir when finished (YES=remove, NO=leave alone)
#   ncp is cp replacement, currently keep as /bin/cp

UNCOMPRESS=gunzip
CLEAN=NO
ncp=/bin/cp


# Given the requested resolution, set dependent resolution parameters
if [[ "$JCAP" = "670" ]]; then
   export LONA=1344
   export LATA=672
   export DELTIM=100
   export resol=1
elif [[ "$JCAP" = "574" ]]; then
   export LONA=1152
   export LATA=576
   export DELTIM=120
   export resol=1
elif [[ "$JCAP" = "382" ]]; then
   export LONA=768
   export LATA=384
   export DELTIM=180
   export resol=1
elif [[ "$JCAP" = "126" ]]; then
   export LONA=384
   export LATA=190
   export DELTIM=600
   export resol=2
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
gdate=`date +%Y%m%d%H -d "${global_4denvar_T126_adate:0:8} ${global_4denvar_T126_adate:8:2} - 6 hours"`
yyg=`echo $gdate | cut -c1-8`
hhg=`echo $gdate | cut -c9-10`
yya=`echo $global_4denvar_T126_adate | cut -c1-8`
hha=`echo $global_4denvar_T126_adate | cut -c9-10`

# Set up $tmpdir
rm -rf $tmpdir
mkdir -p $tmpdir
cd $tmpdir
rm -rf core*

# CO2 namelist and file decisions
ICO2=${ICO2:-0}
if [ $ICO2 -gt 0 ] ; then
        # Copy co2 files to $tmpdir
        co2dir=${CO2DIR:-$fixgsi}
        yyyy=$(echo ${CDATE:-$global_4denvar_T126_adate}|cut -c1-4)
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
        ch4dir=${CH4DIR:-$fixgsi}
        yyyy=$(echo ${CDATE:-$global_4denvar_T126_adate}|cut -c1-4)
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
        n2odir=${N2ODIR:-$fixgsi}
        yyyy=$(echo ${CDATE:-$global_4denvar_T126_adate}|cut -c1-4)
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
        codir=${CODIR:-$fixgsi}
        yyyy=$(echo ${CDATE:-$global_4denvar_T126_adate}|cut -c1-4)
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

# Make gsi namelist

. $scripts/regression_nl_update.sh

SETUP="$SETUP_update"
GRIDOPTS="$GRIDOPTS_update"
BKGVERR="$BKGVERR_update"
ANBKGERR="$ANBKERR_update"
JCOPTS="$JCOPTS_update"
STRONGOPTS="$STRONGOPTS_update"
OBSQC="$OBSQC_update"
OBSINPUT="$OBSINPUT_update"
SUPERRAD="$SUPERRAD_update"
SINGLEOB="$SINGLEOB_update"

if [ "$debug" = ".false." ]; then
   . $scripts/regression_namelists.sh global_4denvar_T126
else
   . $scripts/regression_namelists_db.sh global_4denvar_T126
fi

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
#   satangl  = angle dependent bias correction file (fixed in time)
#   pcpinfo  = text file with information about assimilation of prepcipitation rates
#   ozinfo   = text file with information about assimilation of ozone data
#   errtable = text file with obs error for conventional data (optional)
#   convinfo = text file with information about assimilation of conventional data
#   bufrtable= text file ONLY needed for single obs test (oneobstest=.true.)
#   bftab_sst= bufr table for sst ONLY needed for sst retrieval (retrieval=.true.)

berror=$fixgsi/Big_Endian/global_berror.l${LEVS}y${NLAT}.f77

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
satangl=$fixgsi/global_satangbias.txt
scaninfo=$fixgsi/global_scaninfo.txt
satinfo=$fixgsi/global_satinfo.txt
cloudyinfo=$fixgsi/cloudy_radiance_info.txt
convinfo=$fixgsi/global_convinfo_reg_test.txt
vqcdat=$fixgsi/vqctp001.dat
insituinfo=$fixgsi/global_insituinfo.txt
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

anavinfo=$fixgsi/global_anavinfo.l64.txt
ozinfo=$fixgsi/global_ozinfo.txt
pcpinfo=$fixgsi/global_pcpinfo.txt
errtable=$fixgsi/prepobs_errtable.global
hybens_info=$fixgsi/global_hybens_info.l64.txt
atmsbeamdat=$fixgsi/atms_beamwidth.txt

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
$ncp $scaninfo ./scaninfo
$ncp $satinfo  ./satinfo
$ncp $cloudyinfo  ./cloudy_radiance_info.txt
$ncp $pcpinfo  ./pcpinfo
$ncp $ozinfo   ./ozinfo
$ncp $convinfo ./convinfo
$ncp $vqcdat ./vqctp001.dat
$ncp $insituinfo ./insituinfo
$ncp $errtable ./errtable
$ncp $anavinfo ./anavinfo
$ncp $hybens_info ./hybens_info
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

# Copy CRTM coefficient files based on entries in satinfo file
for file in `awk '{if($1!~"!"){print $1}}' ./satinfo | sort | uniq` ;do
    $ncp $fixcrtm/${file}.SpcCoeff.bin ./
    $ncp $fixcrtm/${file}.TauCoeff.bin ./
done

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

# Copy observational data to $tmpdir
$ncp $global_4denvar_T126_datobs/prepqc.gdas.$global_4denvar_T126_adate                 ./prepbufr
$ncp $global_4denvar_T126_datobs/nsstbufr.gdas.$global_4denvar_T126_adate               ./nsstbufr
$ncp $global_4denvar_T126_datobs/prepbufr.acft_profiles.gdas.$global_4denvar_T126_adate ./prepbufr_profl
$ncp $global_4denvar_T126_datobs/satwnd.gdas.$global_4denvar_T126_adate                 ./satwndbufr
$ncp $global_4denvar_T126_datobs/gpsro.gdas.$global_4denvar_T126_adate                  ./gpsrobufr
$ncp $global_4denvar_T126_datobs/sptrmm.gdas.$global_4denvar_T126_adate                 ./tmirrbufr
$ncp $global_4denvar_T126_datobs/osbuv8.gdas.$global_4denvar_T126_adate                 ./sbuvbufr
$ncp $global_4denvar_T126_datobs/gome.gdas.$global_4denvar_T126_adate                   ./gomebufr
$ncp $global_4denvar_T126_datobs/omi.gdas.$global_4denvar_T126_adate                    ./omibufr
$ncp $global_4denvar_T126_datobs/tcvitl.gdas.$global_4denvar_T126_adate                 ./tcvitl
$ncp $global_4denvar_T126_datobs/goesfv.gdas.$global_4denvar_T126_adate                 ./gsnd1bufr
$ncp $global_4denvar_T126_datobs/1bamua.gdas.$global_4denvar_T126_adate                 ./amsuabufr
$ncp $global_4denvar_T126_datobs/1bamub.gdas.$global_4denvar_T126_adate                 ./amsubbufr
$ncp $global_4denvar_T126_datobs/1bhrs3.gdas.$global_4denvar_T126_adate                 ./hirs3bufr
$ncp $global_4denvar_T126_datobs/1bhrs4.gdas.$global_4denvar_T126_adate                 ./hirs4bufr
$ncp $global_4denvar_T126_datobs/airsev.gdas.$global_4denvar_T126_adate                 ./airsbufr
$ncp $global_4denvar_T126_datobs/mtiasi.gdas.$global_4denvar_T126_adate                 ./iasibufr
$ncp $global_4denvar_T126_datobs/esamua.gdas.$global_4denvar_T126_adate                 ./amsuabufrears
$ncp $global_4denvar_T126_datobs/esamub.gdas.$global_4denvar_T126_adate                 ./amsubbufrears
$ncp $global_4denvar_T126_datobs/eshrs3.gdas.$global_4denvar_T126_adate                 ./hirs3bufrears
$ncp $global_4denvar_T126_datobs/avcsam.gdas.$global_4denvar_T126_adate                 ./avhambufr
$ncp $global_4denvar_T126_datobs/avcspm.gdas.$global_4denvar_T126_adate                 ./avhpmbufr
if [ "$debug" = ".false." ]; then
   $ncp $global_4denvar_T126_datobs/esiasi.gdas.$global_4denvar_T126_adate              ./iasibufrears
fi
$ncp $global_4denvar_T126_datobs/iasidb.gdas.$global_4denvar_T126_adate                 ./iasibufr_db
$ncp $global_4denvar_T126_datobs/gmi1cr.gdas.$global_4denvar_T126_adate                 ./gmibufr
$ncp $global_4denvar_T126_datobs/saphir.gdas.$global_4denvar_T126_adate                 ./saphirbufr
$ncp $global_4denvar_T126_datobs/cris.gdas.$global_4denvar_T126_adate                   ./crisbufr
$ncp $global_4denvar_T126_datobs/crisdb.gdas.$global_4denvar_T126_adate                 ./crisbufr_db
$ncp $global_4denvar_T126_datobs/crisf4.gdas.$global_4denvar_T126_adate                 ./crisfsbufr
$ncp $global_4denvar_T126_datobs/crisf4db.gdas.$global_4denvar_T126_adate               ./crisfsbufr_db
$ncp $global_4denvar_T126_datobs/sevcsr.gdas.$global_4denvar_T126_adate                 ./seviribufr
$ncp $global_4denvar_T126_datobs/atms.gdas.$global_4denvar_T126_adate                   ./atmsbufr
$ncp $global_4denvar_T126_datobs/atmsdb.gdas.$global_4denvar_T126_adate                 ./atmsbufr_db
$ncp $global_4denvar_T126_datobs/ssmisu.gdas.$global_4denvar_T126_adate                 ./ssmisbufr
$ncp $global_4denvar_T126_datobs/abicsr.gdas.$global_4denvar_T126_adate                 ./abibufr
$ncp $global_4denvar_T126_datobs/ahicsr.gdas.$global_4denvar_T126_adate                 ./ahibufr


# Copy bias correction, atmospheric and surface files
$ncp $global_4denvar_T126_datges/biascr.gdas.$gdate                             ./satbias_in
$ncp $global_4denvar_T126_datges/biascr_pc.gdas.${gdate}                        ./satbias_pc
$ncp $global_4denvar_T126_datges/aircraft_t_bias.gdas.$gdate                    ./aircftbias_in
$ncp $global_4denvar_T126_datges/radstat.gdas.$gdate                            ./radstat.gdas

listdiag=`tar xvf radstat.gdas | cut -d' ' -f2 | grep _ges`
for type in $listdiag; do
   diag_file=`echo $type | cut -d',' -f1`
   fname=`echo $diag_file | cut -d'.' -f1`
   date=`echo $diag_file | cut -d'.' -f2`
   $UNCOMPRESS $diag_file
   fnameanl=$(echo $fname|sed 's/_ges//g')
   mv $fname.$date $fnameanl
done

$ncp $global_4denvar_T126_datges/sfnf03.gdas.$gdate  ./sfcf03
$ncp $global_4denvar_T126_datges/sfnf04.gdas.$gdate  ./sfcf04
$ncp $global_4denvar_T126_datges/sfnf05.gdas.$gdate  ./sfcf05
$ncp $global_4denvar_T126_datges/sfnf06.gdas.$gdate  ./sfcf06
$ncp $global_4denvar_T126_datges/sfnf07.gdas.$gdate  ./sfcf07
$ncp $global_4denvar_T126_datges/sfnf08.gdas.$gdate  ./sfcf08
$ncp $global_4denvar_T126_datges/sfnf09.gdas.$gdate  ./sfcf09

$ncp $global_4denvar_T126_datges/nsnf03.gdas.$gdate  ./nstf03
$ncp $global_4denvar_T126_datges/nsnf04.gdas.$gdate  ./nstf04
$ncp $global_4denvar_T126_datges/nsnf05.gdas.$gdate  ./nstf05
$ncp $global_4denvar_T126_datges/nsnf06.gdas.$gdate  ./nstf06
$ncp $global_4denvar_T126_datges/nsnf07.gdas.$gdate  ./nstf07
$ncp $global_4denvar_T126_datges/nsnf08.gdas.$gdate  ./nstf08
$ncp $global_4denvar_T126_datges/nsnf09.gdas.$gdate  ./nstf09

$ncp $global_4denvar_T126_datges/gfngm3.gdas.$global_4denvar_T126_adate  ./sigf03
$ncp $global_4denvar_T126_datges/gfngm2.gdas.$global_4denvar_T126_adate  ./sigf04
$ncp $global_4denvar_T126_datges/gfngm1.gdas.$global_4denvar_T126_adate  ./sigf05
$ncp $global_4denvar_T126_datges/gfnges.gdas.$global_4denvar_T126_adate  ./sigf06
$ncp $global_4denvar_T126_datges/gfngp1.gdas.$global_4denvar_T126_adate  ./sigf07
$ncp $global_4denvar_T126_datges/gfngp2.gdas.$global_4denvar_T126_adate  ./sigf08
$ncp $global_4denvar_T126_datges/gfngp3.gdas.$global_4denvar_T126_adate  ./sigf09

$ncp $global_4denvar_T126_datges/sfcgcy.gdas.$global_4denvar_T126_adate  ./sfcgcy

list="001 002 003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 018 019 020"

for file in $list; do
## ln -s $global_4denvar_T126_datges=/sigf06s_${gdate}_mem${file}_t${JCAP_EN} ./sigf06_ens_mem${file}
   ln -s $global_4denvar_T126_datges/sfg_${gdate}_fhr03s_mem${file} ./sigf03_ens_mem${file}
   ln -s $global_4denvar_T126_datges/sfg_${gdate}_fhr04s_mem${file} ./sigf04_ens_mem${file}
   ln -s $global_4denvar_T126_datges/sfg_${gdate}_fhr05s_mem${file} ./sigf05_ens_mem${file}
   ln -s $global_4denvar_T126_datges/sfg_${gdate}_fhr06s_mem${file} ./sigf06_ens_mem${file}
   ln -s $global_4denvar_T126_datges/sfg_${gdate}_fhr07s_mem${file} ./sigf07_ens_mem${file}
   ln -s $global_4denvar_T126_datges/sfg_${gdate}_fhr08s_mem${file} ./sigf08_ens_mem${file}
   ln -s $global_4denvar_T126_datges/sfg_${gdate}_fhr09s_mem${file} ./sigf09_ens_mem${file}
done

# Run GSI
cd $tmpdir
echo "run gsi now"
eval "$APRUN $tmpdir/gsi.x > stdout 2>&1"
rc=$?
exit $rc
