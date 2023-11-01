!****combine NASA DC-8 merged ICT dataset with modeling data
!   for daily-varied IOAPI files
!
!   Youhua Tang NOAA/ARL and George Mason University

      program conv_group
      include 'PARMS3.EXT'      ! i/o API
      include 'FDESC3.EXT'      ! i/o API
      include 'IODECL3.EXT'     ! i/o API      

      parameter(maxfile=200,nvarmax=900, nfilemax=50)
      
!      parameter(imax=442,jmax=265,kmmax=22,kmax=22,ioff=24,joff=0)
      
      parameter(lchem=71,laero=29,lj=31,laoe=11,ltracer=23)
      
      real, allocatable :: work(:,:,:),topo(:,:),zlevs(:),
     1  zheight(:,:,:),aoe(:,:),dotlon(:,:),dotlat(:,:),conlevs(:),
     2 workdot(:,:,:),workdot2(:,:,:),work2(:,:,:)
 
      real conc(800)
          
!      real work(imax,jmax,kmmax),topo(imax,jmax),zlevs(kmmax),
!     1 zheight(imax,jmax,kmmax),aoe(kmax,laoe),
!     2 dotlon(imax+1,jmax+1),dotlat(imax+1,jmax+1),conlevs(kmmax),
!     3 conc(800),workdot(imax+1,jmax+1,kmmax),projvals(0:8),
!     4 work2(imax,jmax,kmmax)
      
      character afile(20)*80,aline*20000,aline2*400,prefix(2)*80,
     1 suffix(2)*80,name_chem(lchem)*16,name_j(lj)*16,afile2(20)*80,
     2 chtmp*2,afile3(20)*80,name_aero(laero)*16,name_aoe(laoe)*16,
     3 name_tracer(ltracer)*16,gdnam*16,
     4 pcdata(nvarmax)*80
     
      character*80 vname(nvarmax,nfilemax),indexname(nvarmax,nfilemax),
     1 avar(nvarmax),ffile(nfilemax)
      
      integer jday,jtime,nulldate,nulltime,missiondates(23),
     1 julianday(20),btime, nvar(nfilemax),                        ! sampling step in second
     2 indexvar(nvarmax,nfilemax),lcexist(lchem),laexist(laero)
      real pdata(nvarmax),amass(5),convratio(nvarmax)
      logical first,interp4,usefile1,iflag
      real conversion(nvarmax, nfilemax)     ! unit conversion markers
                                             ! >0: conversion factor  =0: no unit conversion
                                             ! -1: ug/m3 to ug/std m3 
                                             ! -2; molecular/cc to ppbv
			                     ! -3; C to K (+273.15)

c      data prefix/'../merge/discoveraq-mrg60-p3b_merge_', 
c     1 'p3b-1m-'/, suffix/'_RB.ict','.dat'/
      data prefix/'../merges/firexaq-mrg60-dc8_merge_', 
     1 'dc8-1m-'/, suffix/'_RL.ict','.dat'/
      data missiondates/20190722,20190724,20190725,20190729,
     1 20190730,20190802,20190803,20190806,20190807,20190808,
     2 20190812,20190813,20190815,20190816,20190819,20190821,
     3 20190823,20190826,20190829,20190830,20190831,20190903,
     4 20190905/
     
      data name_chem/'W_VEL','CO','SO2','SULF','ETHA','O3','NO2','NO',
     1 'PAN','PANX','OPAN','HONO','PNA','NO3','HNO3','NTR1','NTR2',        ! PNA=HNO4, NTR=RNO3
     2 'INTR','N2O5','ETH','FORM','ALD2','ALDX','PAR','OH','HO2','OLE',
     3 'IOLE','H2O2','TOL','XYLMN','ISOP','ISPD','NH3','TERP','FACD',
     4 'AACD','PACD','CRES','CRO','MGLY','MEOH','ETOH','MEPX','ROR',
     5 'C2O3','MEO2','XO2','XO2N','CXO3','HCL','ETHY','CL2','CL',
     6 'ACET','KET','HOCL','CLO','FMCL','BENZENE','RO2',
     7 'BENZRO2','XYLRO2','TOLRO2','OPEN','CRON','CLNO2','CLNO3',
     8 'GLY','GLYD','ROOH'/

      parameter(lnoy=16,lpomi=5,lpomj=7,lsomi=4,lsomj=26)
      character*16 name_noy(lnoy),name_pomi(lpomi),name_pomj(lpomj),
     1 name_somi(lsomi),name_somj(lsomj),
     2 name_tom(lpomi+lsomi+lpomj+lsomj)
      data name_noy/'NO','NO2','NO3','N2O5','HONO','HNO3','PNA','CRON',
     1 'CLNO2','CLNO3','PAN','PANX','OPAN','NTR1','NTR2','INTR'/
      data name_pomi/'ALVPO1I','ASVPO1I','ASVPO2I','APOCI','APNCOMI'/
      data name_pomj/'ALVPO1J','ASVPO1J','ASVPO2J','APOCJ','ASVPO3J',
     1  'AIVPO1J','APNCOMJ'/
      data name_somi/'ALVOO1I','ALVOO2I','ASVOO1I','ASVOO2I'/
      data name_somj/'AISO1J','AISO2J','AISO3J','AMT1J','AMT2J',
     1 'AMT3J','AMT4J','AMT5J','AMT6J','AMTNO3J','AMTHYDJ','AGLYJ',
     2 'ASQTJ','AORGCJ','AOLGBJ','AOLGAJ','ALVOO1J','ALVOO2J','ASVOO1J'
     3 ,'ASVOO2J','ASVOO3J','APCSOJ','AAVB1J','AAVB2J','AAVB3J',
     4 'AAVB4J'/
   
      data name_aero/'ASO4I','ASO4J','ASO4K','ANH4I','ANH4J','ANH4K',
     1 'ANO3I','ANO3J','ANO3K','ANAI','ANAJ','ASEACAT','ACLI','ACLJ',
     2 'ACLK','AECI','AECJ','AOTHRI','AOTHRJ','ACORS','ASOIL','ACAJ',
     3 'AKJ','AFEJ','AALJ','ASIJ','ATIJ','AMNJ','AMGJ'/
      
      data name_j/'NO2_IUPAC10','O3_O3P_IUPAC10','O3_O1D_IUPAC10',
     1 'H2O2_IUPAC10','NO3NO2_06','NO3NO_06','N2O5_IUPAC10',
     2 'HONO_IUPAC10','HNO3_IUPAC10','PNA_IUPAC10','PAN_IUPAC10', 
     3 'MEPX_IUPAC10','NTR_IUPAC10','FORM_R_IUPAC10','FORM_M_IUPAC10', 
     4 'ALD2_R_IUPAC10','ALDX_R_IUPAC10','GLYD_IUPAC10',  
     5 'GLY_R_IUPAC10', 'MGLY_IUPAC10', 'KET_IUPAC10','ACET_IUPAC10', 
     6 'ISPD', 'HPALD','CL2_IUPAC04','HOCL_IUPAC04', 'FMCL_IUPAC04',  
     7  'CLNO2_IUPAC13', 'CLONO2_1', 'CLONO2_2','ACRO_09'/
      
      data name_aoe/'AERO_SSA_W294','EXT_AERO_W294','AERO_SSA_W303',
     1 'EXT_AERO_W303','AERO_SSA_W310','EXT_AERO_W310','AERO_SSA_W381',
     2 'EXT_AERO_W381','AERO_SSA_W607','EXT_AERO_W607','EXT_AERO_W550'/
     
      data name_tracer/'NOX','NOY','SO2','HC1','HC2','HC3','CO',
     1 'SO2LPS','NOxLPS','BIOMCO','CO_TRACER1','CO_TRACER2',
     2 'CO_TRACER3','CO_TRACER4','CO_TRACER5','CO_TRACER6',
     3 'CO_TRACER7','CO_TRACER8','CO_TRACER9','CO_TRACER10',
     4 'CO_TRACER11','CO_TRACER12','CO_TRACER15'/
     
c      data name_cfors/'SO4-AEROSOL','BC','OC','DUST-FINE','DUST-OTHER',
c     1 'CO-Biomass','CO-Fuel','AOD'/ 

      real c303,c302
      parameter(C303=19.83,C302=5417.4)

c --- end WORK_AREAS declarations
      ESAT(TEMK)=.611*EXP(C303-C302/TEMK)       ! for calculating saturated water vapor pressure  
      QSAT(ESAT1,PCB)=ESAT1*.622/(PCB-ESAT1)    ! TEMK is ambient temperature in K, PCB is the pressue in KPa
                                                ! QSAT is the saturated humidity in kg/kg      
      
      name_tom=(/name_pomi,name_somi,name_pomj,name_somj/)
            
      afconst=6.023e23
      paiv=atan(1.)/45 
      radius=6.37e6      ! radius of earth in meter
      
      if(iargc().ne.1) then
        print*,' need flight number'
        stop
       endif
      call getarg(1,aline)
      read(aline,*)nmission 
       
      open(7,file='firex-merge-dc8-rl.ini')   ! variable
      
      nfile=1
      ntotalvar=0
      do while(.true.)
       read(7,'(a800)',end=20)aline
c       print*,trim(aline)
       if(aline(1:1).ne.'#') then   ! comment line
        backspace(7)
        read(7,*)aline,nvar(nfile)
	print*,'nfile,nvar=',nfile,nvar(nfile)
	do i=1,nvar(nfile)
	 read(7,*)vname(i,nfile),indexname(i,nfile),conversion(i,nfile)
c	 print*,vname(i,nfile),indexname(i,nfile),conversion(i,nfile)
	 call lowcase(indexname(i,nfile))
c	 print*,vname(i,nfile),indexname(i,nfile),conversion(i,nfile),
c     1	   i,nfile
	enddo 
	ntotalvar=ntotalvar+nvar(nfile)
	nfile=nfile+1
        endif
       enddo
 20    close(7)
       nfile=nfile-1
       
       print*,'total file and variables are',
     1  nfile,ntotalvar

      if(missiondates(nmission).gt.0) then
       write(aline2,'(i8.8)')missiondates(nmission) 
       open(7,file=trim(prefix(1))//aline2(1:8)//trim(suffix(1)),
     1  status='old')       
      else
       print*,'wrong mission ',nmission,missiondates(nmission)
       stop
      endif
      
      write(aline2,'(i2.2)')nmission
      open(8,file=trim(prefix(2))//aline2(1:2)//trim(suffix(2))) 

101   format('"',a,'",',$)

      do mvar=1,nvar(1)
       write(8,101) trim(vname(mvar,1))
       if(vname(mvar,1).eq.'NO (ppbv)') l_no=mvar
       if(vname(mvar,1).eq.'NO2 (ppbv)') l_no2=mvar
       if(vname(mvar,1).eq.'NOy (ppbv)') l_noy=mvar
       if(vname(mvar,1).eq.'Benzene_WAS (ppbv)') l_benzene=mvar
       if(vname(mvar,1).eq.'Toluene_WAS (ppbv)') l_toluene=mvar
       if(vname(mvar,1).eq.'Latitude (Degs)') l_lat=mvar
       if(vname(mvar,1).eq.'Longitude (Degs)') l_lon=mvar
       if(vname(mvar,1).eq.'Alt_GPS (m)') l_alt=mvar
       if(vname(mvar,1).eq.'Temperature (K)') l_temp=mvar
       if(vname(mvar,1).eq.'RH (%)') l_rh=mvar
       if(vname(mvar,1).eq.'Pressure (Pa)') l_press=mvar
      enddo
      if(l_no.lt.1.or.l_no2.lt.1.or.l_noy.lt.1.or.l_benzene.lt.1.or.
     1 l_toluene.lt.1.or.l_lat.lt.1.or.l_lon.lt.1.or.l_alt.lt.1.or.
     2 l_temp.lt.1.or.l_rh.lt.1.or.l_press.lt.1) then
        print*,'can not find index' 
	print*,l_no,l_no2,l_noy,l_benzene,l_toluene,l_lat,l_lon,
     1    l_alt,l_temp,l_rh,l_press	
	stop
       endif
      
      read(7,*)nline,itmp ! header line number
      nline=nline-1
      do i=1,5                  ! skip 5 line
       read(7,*)
      enddo
         
      read(7,*) iyear,imonth,idate                          ! date information is the 7st line
      jday=iyear*1000+julian(iyear,imonth,idate)            ! YYYYDDD

      nline=nline-7
         
      do i=1,nline
        read(7,*)
      enddo
      read(7,'(a20000)')aline
      call decommach(aline,avar,numavarfile)

      do n=1,nvar(1)
       do i=1,numavarfile
        call lowcase(avar(i))
	if(trim(indexname(n,1)).eq.trim(adjustl(avar(i)))) exit
       enddo
       indexvar(n,1)=i 	
       
       if(i.gt.numavarfile) then
        print*,'can not find ',trim(vname(n,1))
	print*,trim(aline)
	do i=1,numavarfile
	  print*,'i,numavarfile,avar(i)=',
     1	    i,numavarfile,avar(i)
        enddo
	indexvar(n,1)=-999
       endif
      enddo  	
	 	 
      write(8,102)
 102  format('"NOx (ppbv)","NOz (ppbv)","Benzene+Toluene (ppbv)",',
     2 '"model_xloc","model_yloc",',
     3 '"model_WS (m/s)","model_WD","model_Temp (K)","model_RH (%)",',  ! meteo
     4 '"model_P (Pa)","model_Vapor (g/kg)","model_CWATER (g/kg)",',
     4 '"model_W (m/s)","model_CO (ppbv)",',         ! chem
     5 '"model_SO2 (ppbv)","model_H2SO4 (ppbv)","model_C2H6 (ppbv)",',
     6 '"model_O3 (ppbv)","model_NO2 (ppbv)","model_NO (ppbv)",',
     7 '"model_NOx (ppbv)","model_PAN (ppbv)","model_PANX (ppbv)",',
     7 '"model_OPAN (ppbv)",',
     7 '"model_HONO (ppbv)","model_HNO4 (ppbv)","model_NO3 (ppbv)",',
     8 '"model_HNO3 (ppbv)","model_NTR1 (ppbv)","model_NTR2 (ppbv)",',
     8 '"model_INTR (ppbv)","model_N2O5 (ppbv)",',
     9 '"model_C2H4 (ppbv)",',
     a '"model_HCHO (ppbv)","model_CCHO (ppbv)","model_ALDX (ppbv)",',
     b '"model_PAR (ppbv)","model_OH (ppbv)","model_HO2 (ppbv)",',
     c '"model_OLE (ppbv)","model_IOLE (ppbv)",',
     c '"model_H2O2 (ppbv)","model_TOL (ppbv)","model_XYLMN (ppbv)",',
     d '"model_ISOP (ppbv)","model_ISPD (ppbv)","model_NH3 (ppbv)",',
     e '"model_TERP (ppbv)","model_FACD (ppbv)","model_AACD (ppbv)",',
     f '"model_PACD (ppbv)","model_CRES (ppbv)","model_CRO (ppbv)",',
     g '"model_MGLY (ppbv)","model_MEOH (ppbv)","model_ETOH (ppbv)",',
     h '"model_MEPX (ppbv)","model_ROR (ppbv)","model_C2O3 (ppbv)",',
     i '"model_MEO2 (ppbv)","model_XO2 (ppbv)","model_XO2N (ppbv)",',
     j '"model_CXO3 (ppbv)","model_HCL (ppbv)","model_ETHY (ppbv)",',
     k '"model_CL2 (ppbv)","model_CL (ppbv)","model_ACET (ppbv)",',
     l '"model_KET (ppbv)","model_HOCL (ppbv)","model_CLO (ppbv)",',
     m '"model_FMCL (ppbv)","model_BENZENE (ppbv)","model_RO2 (ppbv)"',
     n ',"model_BENZRO2 (ppbv)","model_XYLRO2 (ppbv)",',
     o '"model_TOLRO2 (ppbv)","model_OPEN (ppbv)","model_CRON (ppbv)"',
     o ',"model_CLNO2 (ppbv)","model_CLNO3 (ppbv)",',
     p '"model_NOy (ppbv)","model_NOz (ppbv)","model_GLY (ppbv)",',
     p '"model_GLYD (ppbv)","model_ROOH (ppbv)",', 
     p '"model_PM1AT","model_PM1AC",',
     p '"model_SO4I (ug/std m3)","model_SO4J (ug/std m3)",',          ! Sulfate
     q '"model_SO4K (ug/std m3)","model_PM1_SO4 (ug/std m3)",',
     q '"model_tot_SO4 (ug/std m3)","SO4_Fine_ratio",',
     p '"model_NH4I (ug/std m3)","model_NH4J (ug/std m3)",',          ! NH4
     q '"model_NH4K (ug/std m3)","model_PM1_NH4 (ug/std m3)",',
     q '"model_tot_NH4 (ug/std m3)","NH4_Fine_ratio",',
     p '"model_NO3I (ug/std m3)","model_NO3J (ug/std m3)",',          ! Nitrate
     q '"model_NO3K (ug/std m3)","model_PM1_NO3 (ug/std m3)",',
     q '"model_tot_NO3 (ug/std m3)","NO3_Fine_ratio",',
     p '"model_NaI (ug/std m3)","model_NaJ (ug/std m3)",',          ! sodium
     q '"model_NaK (ug/std m3)","model_PM1_Na (ug/std m3)",',
     q '"model_tot_Na (ug/std m3)","Na_Fine_ratio",',
     p '"model_ClI (ug/std m3)","model_ClJ (ug/std m3)",',          ! chloride
     q '"model_ClK (ug/std m3)","model_PM1_Cl (ug/std m3)",',
     q '"model_tot_Cl (ug/std m3)","Cl_Fine_ratio",',
     r '"model_ECI (ug/std m3)","model_ECJ (ug/std m3)",',  ! EC
     r '"model_PM1_EC (ug/std m3)",', 
     r '"model_OTHRI (ug/std m3)","model_OTHRJ (ug/std m3)",', ! Other
     r '"model_PM1_OTHR (ug/std m3)",',
     s '"model_ACORS (ug/std m3)","model_ASOIL (ug/std m3)",',
     s '"model_CAJ (ug/std m3)","model_PM1_CA (ug/std m3)",',
     t '"model_KJ (ug/std m3)","model_PM1_K (ug/std m3)",',
     t '"model_FEJ (ug/std m3)","model_PM1_FE (ug/std m3)",',
     t '"model_ALJ (ug/std m3)","model_PM1_AL (ug/std m3)",',
     u '"model_SIJ (ug/std m3)","model_PM1_SI (ug/std m3)",',
     v '"model_TIJ (ug/std m3)","model_PM1_TI (ug/std m3)",',
     v '"model_MNJ (ug/std m3)","model_PM1_MN (ug/std m3)",',
     x '"model_MGJ (ug/std m3)","model_PM1_MG (ug/std m3)",',
     p '"model_tot_NO3 (gas+aerosol) (ppbv)",',     
     y '"model_POAI (ug/std m3)","model_SOAI (ug/std m3)",',   ! organic aerosols
     y '"model_TOAI (ug/std m3)","model_POAJ (ug/std m3)",',
     z '"model_SOAJ (ug/std m3)","model_TOAJ (ug/std m3)",',
     z '"model_PM1_POA (ug/std m3)","model_PM1_SOA (ug/std m3)",',
     z '"model_PM1_TOA (ug/std m3)",'
     1 '"model_J-NO2 (/s)","model_J-O3P (/s)","model_J-O1D (/s)",', ! jvalues
     2 '"model_J-H2O2 (/s)","model_J-NO3-NO2 (/s)",',
     3 '"model_J-NO3-NO (/s)","model_J-N2O5 (/s)","model_J-HONO (/s)"',
     4 ',"model_J-HNO3 (/s)","model_J-HNO4 (/s)","model_J-PAN (/s)",',
     5 '"model_J-ROOH (/s)","model_J-NTR (/s)","model_J-HCHO-R (/s)",',
     6 '"model_J-HCHO-M (/s)","model_J-ALD2 (/s)","model_J-ALDX (/s)"',
     7 ',"model_J-GLYD (/s)","model_J-GLY (/s)","model_J-MGLY (/s)",',
     8 '"model_J-KET (/s)","model_J-ACET (/s)","model_J-ISPD (/s)",',
     9 '"model_J-HPALD (/s)","model_J-CL2 (/s)","model_J-HOCL (/s)",',
     a '"model_J-FMCL (/s)","model_J-CLNO2 (/s)",',
     b '"model_J-CLNO3-NO2 (/s)","model_J-CLNO3-NO3 (/s)",',
     c '"model_J-ACRO (/s)",',
     d '"model_SSA@294","model_AOE@294 (km-1)","model_SSA@303",',
     e '"model_AOE@303 (km-1)","model_SSA@310","model_AOE@310 (km-1)"',
     f ',"model_SSA@381","model_AOE@381 (km-1)","model_SSA@607",',
     g '"model_AOE@607 (km-1)","model_AOE@550 (km-1)"')     ! extinction 
     
      if(.not.open3('TOPO',FSREAD3,'pathway')) then   
       print*,' Error open TOPO file'
       stop
      endif
      if(.not.desc3('TOPO')) STOP
      imax=ncols3d
      jmax=nrows3d
      
      if(index(gdnam3d,'CON_12').ge.1.or.
     1   index(gdnam3d,'CONUS').ge.1) then
        gdnam='AQF_CONUS_5x'
      else if (index(gdnam3d,'COL_04').ge.1) then
        gdnam='COL_04'
      else
       print*,'unknown grid ',gdnam3d
       stop
      endif

      ddx=sngl(xcell3d)
      xorig=sngl(xorig3d)
      yorig=sngl(yorig3d)
       		
      allocate(topo(imax,jmax))            
      if(.not.READ3('TOPO','HT',ALLAYS3, 0, 0, topo)) then
        print*, 'Error in reading topograph'
	stop
      endif      
      iflag=close3('TOPO')
c----dot files
       allocate(dotlon(imax+1,jmax+1),dotlat(imax+1,jmax+1))
      if(.not.open3('TOPODOT',FSREAD3,'pathway')) then
       print*,' Error open TOPODOT file'
       stop
      endif
      if(.not.READ3('TOPODOT','LATD',ALLAYS3, 0, 0, dotlat)) then
        print*, 'Error in reading LAT'
        stop
      endif
      if(.not.READ3('TOPODOT','LOND',ALLAYS3, 0, 0, dotlon)) then
        print*, 'Error in reading LON'
        stop
      endif	
      iflag=close3('TOPODOT')

      if(.not.open3('WINDDOT',FSREAD3,'stem2v5d')) stop
      if (.not.DESC3('WINDDOT') ) stop   ! get grid information
      
      if(ncols3d.ne.imax+1.or.nrows3d.ne.jmax+1) then
       print*,'dimension wrong ',ncols3d,nrows3d
       stop
      endif
      if(.not.open3('WINDDOT2',FSREAD3,'stem2v5d')) stop

c-----cross point files


      if(.not.open3('MET3DCRO',FSREAD3,'stem2v5d')) then
       print*,' Error open MET3DCRO file'
       stop
      endif

      if(.not.open3('MET3DCRO2',FSREAD3,'stem2v5d')) stop

      if (.not. DESC3('MET3DCRO') ) stop   ! get grid information from 3d chemical output
      
      kmmax=nlays3d
      allocate(zheight(imax,jmax,kmmax),work(imax,jmax,kmmax))
      allocate(work2(imax,jmax,kmmax),workdot(imax+1,jmax+1,kmmax),
     1 workdot2(imax+1,jmax+1,kmmax) )
      allocate(conlevs(kmmax),zlevs(kmmax))
      
      if(.not.READ3('MET3DCRO','ZH',ALLAYS3,sdate3d,stime3d,zheight))
     1  then
        print*, 'Error in reading MET3DCRO'
        stop
      endif

      if(ncols3d.ne.imax.or.nrows3d.ne.jmax) then
       print*,'dimension wrong ',ncols3d,nrows3d
       stop
      endif

      do i=1,imax
       do j=1,jmax
        do k=1,kmmax
        zheight(i,j,k)=topo(i,j)+zheight(i,j,k)    ! vertical altitude
        enddo
       enddo
      enddo

      if(.not.OPEN3('CHEM3D',FSREAD3,'pathway')) then
       print*,'open input file error for CHEM3D'
       stop
      endif

      if (.not. DESC3('CHEM3D') ) then   ! get grid information from 3d chemical
       print*, 'Error getting info from CHEM3D'
       stop
      endif

      lcexist(1:lchem)=0
      laexist(1:laero)=0
      do j=1,nvars3d
       do i=1,lchem
        if(name_chem(i).eq.vname3d(j)) lcexist(i)=1
       enddo
       do i=1,laero
        if(name_aero(i).eq.vname3d(j)) laexist(i)=1
       enddo	
      enddo
        
      if(ncols3d.ne.imax.or.nrows3d.ne.jmax.or.nlays3d.gt.kmmax) then
       print*,'dimension wrong ',ncols3d,nrows3d,nlays3d
       stop
      endif

      kmax=nlays3d
      kchem3d=nlays3d
      ichemstep=tstep3d/10000
       
      ichem1_hour_end=(sdate3d-iyear*1000)*24+
     1  (stime3d+tstep3d*(mxrec3d-1))/10000            ! end julian hour in CHEM1
     
      print*,'ichem1_hour_end,sdate3d,stime3d,tstep3d,mxrec3d=',
     1  ichem1_hour_end,sdate3d,stime3d,tstep3d,mxrec3d

      if(.not.OPEN3('CHEM3D2',FSREAD3,'pathway')) then
	print*,'open input file error for CHEM3D2'
       stop
      endif

      if(.not.OPEN3('PMDIAG',FSREAD3,'pathway')) stop
      if(.not.OPEN3('PMDIAG2',FSREAD3,'pathway')) stop
      
      if(.not.OPEN3('JVFILE',FSREAD3,'pathway')) stop
      if(.not.OPEN3('JVFILE2',FSREAD3,'pathway')) stop
      
      if(.not.OPEN3('AOP',FSREAD3,'pathway')) stop
      if(.not.OPEN3('AOP2',FSREAD3,'pathway')) stop
      
      IF(.NOT.LAMBERT(GDNAM, P_ALP3d, P_BET3d, P_GAM3d,
     1       XCENT3D, YCENT3D )) stop

      do while(.true.) 
       read(7,'(a20000)',end=99)aline
!       print*,'read line ',trim(aline)
       call decomma(aline,pdata,ndata)
       if(ndata.ne.numavarfile) then
        print*,ndata,numavarfile,'inconsistent varnumber'
	print*,trim(aline)
	stop
       endif
       
       do L=1,nvar(1)
         i=indexvar(L,1)
	 
	 if(i.gt.0) then
           if(conversion(L,1).gt.1e-19.and.
     1	       pdata(i).gt.-888) then
	     conc(L)=pdata(i)*conversion(L,1)
 	   elseif(abs(conversion(L,1)+1).le.1e-19.and.  ! conversion = -1: ug/m3 to ug/std m3 
     1	     pdata(i).gt.-888) then
            conc(L)=pdata(i)*temp/press*101300/273 
           elseif(abs(conversion(L,1)+2).le.1e-19.and.      ! conversion = -2: molecular/cc to ppbv
     1	          pdata(i).gt.-888) then	
            conc(L)=pdata(i)*temp/press/7.2427e7
           elseif(abs(conversion(L,1)+3).le.1e-19.and.      ! conversion = -3: c to K
     1	          pdata(i).gt.-888) then	
            conc(L)=pdata(i)+273.15
 	   elseif(abs(conversion(L,1)+4).le.1e-19.and.  ! conversion = -4: ng/m3 to ug/std m3 
     1	     pdata(i).gt.-888) then
            conc(L)=pdata(i)*temp/press*101300/273./1000.
 	   elseif(abs(conversion(L,1)+5).le.1e-19.and.  ! conversion = -4: ng/m3 to ug/std m3 
     1	     pdata(i).gt.-888) then
            pdata(i)=pdata(i)+30
            conc(L)=pdata(i)/3600.     ! Time+30second -> hour
           elseif(abs(conversion(L,1)).le.1e-19) then    ! no conversion
	     conc(L)=pdata(i)
           else
             if(pdata(i).gt.-888) then
               print*,'wrong conversion ',vname(L,1),
     1		 conversion(L,1)
               stop
	     else
		conc(L)=pdata(i)
	     endif 
	   endif
	   if(L.eq.l_temp) temp=conc(L)
	   if(L.eq.l_press) press=conc(L)
	 else
	  conc(L)=-9999.
	 endif         
       	 
	if(vname(L,1).eq.'JDAY') then
	 jday=int(pdata(i))
	 write(8,"(i3,$)")jday
	else if(vname(L,1).eq.'UTC (hour)') then
	 write(8,"(f9.6,$)")conc(L) 
	 ihour=int(pdata(i))/3600                ! pdata(i) is in second
	 iminute=(int(pdata(i))-ihour*3600)/60
	 isecond=mod(int(pdata(i)),60)
	else  
	 if(conc(L).gt.-888) write(8,"(g12.5,$)")conc(L)
	endif 
        write(8,"(',',$)")
       enddo 
	         
       cnoy=0.
       cnoz=0.
       Lmodeltotal=6+lchem+3+laero+3+4*3                       ! meteorology, chem, NOy+NOz+NOx, aerosol, total BC+A25+NH4
     	                                                      ! total_ions+F_ion+ion_fine ratios, total_NO3

       hourjulian=jday*24+ihour+iminute/60.+isecond/3600.    ! julian hour
               
       if(ihour.ge.24) then
        ihour=ihour-24
	jday=jday+1
       endif	
       
       nowtime=ihour*10000+iminute*100+isecond   ! HHMMSS
       
       print*,'ichem1_hour_end,hourjulian,jday,nowtime=',
     1 ichem1_hour_end,hourjulian,jday,nowtime              

       
       if(conc(l_no).ge.0.and.conc(l_no2).ge.0) then     ! write observed NOx
         vnox=conc(l_no)+conc(l_no2)
	 write(8,'(g10.5,$)')vnox
       else
         vnox=-9999.
       endif
       
       write(8,"(',',$)")
       
       if(vnox.ge.0.and.conc(l_noy).ge.0)            ! write NOz
     1     write(8,'(g10.5,$)') conc(l_noy)-vnox

       write(8,"(',',$)")       
       if(conc(l_benzene).gt.0.and.conc(l_toluene).gt.0)
     1  write(8,'(g10.5,$)')conc(l_benzene)+conc(l_toluene)
       
***Begin interpolating model value

       if(hourjulian.gt.ichem1_hour_end) then
        usefile1=.false.
       else
        usefile1=.true.
       endif  
       
       if(conc(l_lon).le.-900.or.conc(l_lat).le.-900.or.
     1     conc(l_alt).le.-900) then
         print*,'no valid lat/lon/alt there '
         write(8,'(1x)')
         cycle
        endif

        if(.not.LL2LAM(conc(l_lon),conc(l_lat),x,y)) stop
	x=(x-xorig)/ddx+0.5
	y=(y-yorig)/ddx+0.5

	write(8,"(2(',',f9.3),$)")x,y     ! print out grid locations related to domain
	
	if(x.gt.imax.or.x.lt.1.or.y.gt.jmax.or.y.lt.1.) then
	 do L=1,Lmodeltotal
	  write(8,"(',',$)")
	 enddo
	 print*,'out of model domain ',conc(l_lat),conc(l_lon),x,y,
     1	  'skip ',Lmodeltotal,' species'
         write(8,'(1x)')
	 cycle
	endif

	xratio=x-int(x)
	yratio=y-int(y)
        
	xd=x-0.5
        yd=y-0.5
        xdratio=xd-int(xd)
        ydratio=yd-int(yd)
C---zlevs

      do mlev=1,kmax
	zlevs(mlev)=(1-yratio)*(zheight(int(x),int(y),mlev)*
     1	 (1-xratio)+zheight(int(x)+1,int(y),mlev)*xratio)+yratio*
     2   (zheight(int(x),int(y)+1,mlev)*(1-xratio)+
     3    zheight(int(x)+1,int(y)+1,mlev)*xratio)
       enddo

       fheight=conc(l_alt)
       do mlev=1,kmax
        if(fheight.le.zlevs(mlev)) exit
       enddo
       kp=mlev
       if(kp.eq.1) then    !if below the lowest layer
       	kp=2
        zratio=0.
       else if(fheight.gt.zlevs(mlev)) then !if above top
        zratio=1.
	kp=kmax
       else
        zratio=(fheight-zlevs(kp-1))/
     1   (zlevs(kp)-zlevs(kp-1))
       endif

	print*,'lat,lon,x,y, fheight,kp, zratio, nowtime=',conc(l_lat),
     1	conc(l_lon),x,y,fheight,kp,zratio, nowtime, usefile1
        print*,'imax,jmax,kmax,kmmax,iyear,jday=',imax,jmax,kmax,kmmax,
     1	 iyear,jday
	
!-------U
       if(usefile1) then
        if(.not.INTERP3('WINDDOT','UWIND','pathway',
     1   iyear*1000+jday,nowtime,(imax+1)*(jmax+1)*kmax,workdot)) stop 
       else  
        if(.not.INTERP3('WINDDOT2','UWIND','pathway',
     1   iyear*1000+jday,nowtime,(imax+1)*(jmax+1)*kmax,workdot)) stop 
       endif
       print*,'horizontal interpolate U'
       	 
       do mlev=1,kmax
	conlevs(mlev)=(1-ydratio)*(workdot(int(xd),int(yd),mlev)*
     1	(1-xdratio)+workdot(int(xd)+1,int(yd),mlev)*xdratio)+
     2  ydratio*(workdot(int(xd),int(yd)+1,mlev)*(1-xdratio)+
     3   workdot(int(xd)+1,int(yd)+1,mlev)*xdratio)
       end do 
	u=conlevs(kp)*zratio+conlevs(kp-1)*(1-zratio)
   	 
C** START reading modeling meteorology

       if(usefile1) then
        if(.not.INTERP3('WINDDOT','VWIND','pathway',
     1   iyear*1000+jday,nowtime,(imax+1)*(jmax+1)*kmax,workdot2)) stop 
       else  
        if(.not.INTERP3('WINDDOT2','VWIND','pathway',
     1   iyear*1000+jday,nowtime,(imax+1)*(jmax+1)*kmax,workdot2)) stop 
       endif
       
	do mlev=1,kmax
	 conlevs(mlev)=(1-ydratio)*(workdot2(int(xd),int(yd),mlev)*
     1	 (1-xdratio)+workdot2(int(xd)+1,int(yd),mlev)*xdratio)+ydratio*
     2   (workdot2(int(xd),int(yd)+1,mlev)*(1-xdratio)+
     3    workdot2(int(xd)+1,int(yd)+1,mlev)*xdratio)
        end do 
	 v=conlevs(kp)*zratio+conlevs(kp-1)*(1-zratio)
	 ws=sqrt(u*u+v*v)         ! modeling wind speed
	 
	 write(8,"(',',f5.2,$)")ws
	 if(u.eq.0) then   ! calculate wind direction in MET coordinate  
	  if(v.lt.0) then
 	   rams_wd=0.
	  else
	   rams_wd=180. 
	  endif 
	 else 
	  rams_wd=atan(v/u)/paiv
	  if(u.lt.0) then
	   rams_wd=90.-rams_wd
	  else
	   rams_wd=270.-rams_wd 
	  endif
	 endif              
	 wd=rams_wd-asin(radius*cos(dotlat(nint(xd),nint(yd+1))*paiv)        ! wd in geographical coordinate
     1    *paiv*(dotlon(nint(xd),nint(yd+1))-dotlon(nint(xd),nint(yd)))
     2    /ddx)/paiv
	 if(wd.lt.0) wd=wd+360.  
	 write(8,"(',',f6.2,$)")wd
        
	print*,'read met3dcro'
	
	if(usefile1) then
 	 if(.not.INTERP3('MET3DCRO','TA','meteo',             ! reading temperature
     1	  iyear*1000+jday,nowtime,imax*jmax*kmmax,work)) stop
        else
 	 if(.not.INTERP3('MET3DCRO2','TA','meteo',
     1	  iyear*1000+jday,nowtime,imax*jmax*kmmax,work)) stop	
	endif 
	
	 do mlev=1,kmax
	 conlevs(mlev)=(1-yratio)*(work(int(x),int(y),mlev)*
     1	 (1-xratio)+work(int(x)+1,int(y),mlev)*xratio)+yratio*
     2   (work(int(x),int(y)+1,mlev)*(1-xratio)+
     3    work(int(x)+1,int(y)+1,mlev)*xratio)
         end do 
	 temp=conlevs(kp)*zratio+conlevs(kp-1)*(1-zratio)

	if(usefile1) then
 	 if(.not.INTERP3('MET3DCRO','PRES','meteo',             ! reading pressure
     1	  iyear*1000+jday,nowtime,imax*jmax*kmmax,work)) stop
        else
 	 if(.not.INTERP3('MET3DCRO2','PRES','meteo',
     1	  iyear*1000+jday,nowtime,imax*jmax*kmmax,work)) stop	
	endif 
	 do mlev=1,kmax
	 conlevs(mlev)=(1-yratio)*(work(int(x),int(y),mlev)*
     1	 (1-xratio)+work(int(x)+1,int(y),mlev)*xratio)+yratio*
     2   (work(int(x),int(y)+1,mlev)*(1-xratio)+
     3    work(int(x)+1,int(y)+1,mlev)*xratio)
         end do 
	 press=conlevs(kp)*zratio+conlevs(kp-1)*(1-zratio)

	if(usefile1) then
 	 if(.not.INTERP3('MET3DCRO','QV','meteo',             ! water vapor
     1	  iyear*1000+jday,nowtime,imax*jmax*kmmax,work)) stop
        else
 	 if(.not.INTERP3('MET3DCRO2','QV','meteo',
     1	  iyear*1000+jday,nowtime,imax*jmax*kmmax,work)) stop	
	endif 
	 do mlev=1,kmax
	 conlevs(mlev)=(1-yratio)*(work(int(x),int(y),mlev)*
     1	 (1-xratio)+work(int(x)+1,int(y),mlev)*xratio)+yratio*
     2   (work(int(x),int(y)+1,mlev)*(1-xratio)+
     3    work(int(x)+1,int(y)+1,mlev)*xratio)
         end do 
	 vapor=conlevs(kp)*zratio+conlevs(kp-1)*(1-zratio)
	 rh=100.*vapor/QSAT(ESAT(         ! computing relative humudity
     1             temp),press/1000)           !convert to KPa

         write(8,"(2(',',f6.2),',',f9.1,',',g10.5,$)")temp,rh,press,
     1	   vapor*1000 

	if(usefile1) then
 	 if(.not.INTERP3('MET3DCRO','QC','meteo',             ! cloud water
     1	  iyear*1000+jday,nowtime,imax*jmax*kmmax,work)) stop
        else
 	 if(.not.INTERP3('MET3DCRO2','QC','meteo',
     1	  iyear*1000+jday,nowtime,imax*jmax*kmmax,work)) stop	
	endif 
	 do mlev=1,kmax
	 conlevs(mlev)=(1-yratio)*(work(int(x),int(y),mlev)*
     1	 (1-xratio)+work(int(x)+1,int(y),mlev)*xratio)+yratio*
     2   (work(int(x),int(y)+1,mlev)*(1-xratio)+
     3    work(int(x)+1,int(y)+1,mlev)*xratio)
         end do 
         cwater= conlevs(kp)*zratio+conlevs(kp-1)*(1-zratio)  ! Cloud Water content in Kg/Kg

	if(usefile1) then
 	 if(.not.INTERP3('MET3DCRO','QR','meteo',             ! Rain water
     1	  iyear*1000+jday,nowtime,imax*jmax*kmmax,work)) stop
        else
 	 if(.not.INTERP3('MET3DCRO2','QR','meteo',
     1	  iyear*1000+jday,nowtime,imax*jmax*kmmax,work)) stop	
	endif 
	 do mlev=1,kmax
	 conlevs(mlev)=(1-yratio)*(work(int(x),int(y),mlev)*
     1	 (1-xratio)+work(int(x)+1,int(y),mlev)*xratio)+yratio*
     2   (work(int(x),int(y)+1,mlev)*(1-xratio)+
     3    work(int(x)+1,int(y)+1,mlev)*xratio)
         end do 
         rwater=conlevs(kp)*zratio+conlevs(kp-1)*(1-zratio)  ! Rain Water content in Kg/Kg
!	 cwater=1000*(cwater+rwater)/1.28                   ! air density STD in kg/m3
                                                           ! convert to g/std m3
!	 write(8,"(',',g10.5,$)")cwater*press/temp*273./101300  ! g/std m3 to g/m3
         write(8,"(',',g10.5,$)")(cwater+rwater)*1000
	 
C****Starting Modeling Gas-phase
         
	 do L=1,lchem

        if(lcexist(L).eq.1) then
         if(hourjulian.le.ichem1_hour_end) then	  	 
 	  if(L.le.1) print*,'read CHEM3D ', hourjulian,ichem1_hour_end
	  if(.not.INTERP3('CHEM3D',name_chem(L),'pathway',
     1	   iyear*1000+jday,nowtime,imax*jmax*kchem3d,work)) then
           print*,'interp error for ',name_chem(L),' on ',jday,nowtime	 

c   if(name_chem(L).eq.'CO2'.or.name_chem(L).eq.'METHACRO'
c    1	   .or. name_chem(L).eq.'MEOH'.or.name_chem(L).eq.'HCOOH'
c    2     .or. name_chem(L).eq.'CCO_OH'.or.name_chem(L).eq.'RCHO') then
c    write(8,"(',',$)")           ! skip these species
c    goto 56
c   endif 
	   stop
	  endif 
	  
!	 elseif(hourjulian.ge.(ichem1_hour_end+ichemstep)) then 
         else
	  if(L.le.1) print*,'read CHEM3D2', hourjulian,ichem1_hour_end
	  if(.not.INTERP3('CHEM3D2',name_chem(L),'pathway',
     1	   iyear*1000+jday,nowtime,imax*jmax*kchem3d,work)) then
           print*,'interp error for ',name_chem(L),' on ',jday,nowtime	 
c   if(name_chem(L).eq.'CO2'.or.name_chem(L).eq.'METHACRO'
c    1	   .or. name_chem(L).eq.'MEOH'.or.name_chem(L).eq.'HCOOH'
c    2     .or. name_chem(L).eq.'CCO_OH'.or.name_chem(L).eq.'RCHO') then
c    write(8,"(',',$)")           ! skip these species
c    goto 56
c   endif 
	   stop
	  endif 
	 
!	 else          !    ichem1_hour_end  <  hourjulian< (ichem1_hour_end+ichemstep)
!	  idaytmp1=ichem1_hour_end/24
!	  itimetmp1=mod(ichem1_hour_end,24)
	  
!	  if(.not.READ3('CHEM3D',name_chem(L),ALLAYS3,
!     1      iyear*1000+idaytmp1,itimetmp1*10000,work)) stop

!	  idaytmp2=(ichem1_hour_end+ichemstep)/24
!	  itimetmp2=mod(ichem1_hour_end+ichemstep,24)	 
!	  if(.not.READ3('CHEM3D2',name_chem(L),ALLAYS3,
!     1      iyear*1000+idaytmp2,itimetmp2*10000,work2)) stop      ! read the time step of the second file
	  
!	  tratio=(hourjulian-ichem1_hour_end)/ichemstep
!	  if(tratio.lt.0.or.tratio.gt.1) then
!	   print*,'wrong tratio',hourjulian,ichem1_hour_end,ichemstep
!	   stop
!	  endif 
!	  work(1:imax,1:jmax,1:kmax)=work(1:imax,1:jmax,1:kmax)
!     1     *(1-tratio) +  work2(1:imax,1:jmax,1:kmax)*tratio
	  
	 endif
	 
	 do mlev=1,kmax
	  conlevs(mlev)=(1-yratio)*(work(int(x),int(y),mlev)*
     1	  (1-xratio)+work(int(x)+1,int(y),mlev)*xratio)+yratio*
     2    (work(int(x),int(y)+1,mlev)*(1-xratio)+
     3     work(int(x)+1,int(y)+1,mlev)*xratio)
         enddo

        conc(L)=(conlevs(kp)*zratio+conlevs(kp-1)*(1-zratio))
	if(name_chem(L).ne.'W_VEL') conc(L)=conc(L)*1000.   ! ppmV to ppbV
	
        if(conc(L).le.0.01) then
	  write(8,"(',',g13.5,$)")conc(L) 
	else
	  write(8,"(',',g10.5,$)")conc(L)
	endif   
	
       else  ! lcexist=0
        write(8,"(',',$)")
       endif	
	
	if(name_chem(L).eq.'NO') 
     1	  write(8,"(',',g10.5,$)")conc(L)+conc(L-1) !NOx

        if(name_chem(L).eq.'HNO3') hno3=conc(L)   ! store HNO3
	if(any(name_noy.eq.name_chem(L))) then
           cnoy=cnoy+conc(L)                     ! accumulate NOy
	   if(name_chem(L).ne.'NO'.and.name_chem(L).ne.'NO2')
     1      cnoz=cnoz+conc(L)                     ! accumulate NOz
	   if(name_chem(L).eq.'N2O5') then 
	    cnoy=cnoy+conc(L)
	    cnoz=cnoz+conc(L)
	  endif	  
        endif 
 56    if(name_chem(L).eq.'CLNO3') then
	 if(cnoy.ge.0) then
 	  write(8,"(2(',',f8.3),$)")cnoy,cnoz ! write out NOy, NOz after CLNO3
	 else
	  write(8,"(',,',$)")
	 endif  
	 cnoy=0.
	 cnoz=0.
	endif 

        enddo

C** reading PM1 fraction
	if(hourjulian.le.ichem1_hour_end) then
	 aline2='PMDIAG'
	else
	 aline2='PMDIAG2'
        endif	 
	if(.not.interp3(trim(aline2),'PM1AT','pathway',
     1	  iyear*1000+jday,nowtime,imax*jmax*kchem3d,work)) stop
      
     	 do mlev=1,kchem3d
	  conlevs(mlev)=(1-yratio)*(work(int(x),int(y),mlev)*
     1	  (1-xratio)+work(int(x)+1,int(y),mlev)*xratio)+yratio*
     2    (work(int(x),int(y)+1,mlev)*(1-xratio)+
     3     work(int(x)+1,int(y)+1,mlev)*xratio)
         enddo
         pm1at=conlevs(kp)*zratio+conlevs(kp-1)*(1-zratio)

	 if(.not.interp3(trim(aline2),'PM1AC','pathway',
     1	  iyear*1000+jday,nowtime,imax*jmax*kchem3d,work)) stop 
     	 do mlev=1,kchem3d
	  conlevs(mlev)=(1-yratio)*(work(int(x),int(y),mlev)*
     1	  (1-xratio)+work(int(x)+1,int(y),mlev)*xratio)+yratio*
     2    (work(int(x),int(y)+1,mlev)*(1-xratio)+
     3     work(int(x)+1,int(y)+1,mlev)*xratio)
         enddo
         pm1ac=conlevs(kp)*zratio+conlevs(kp-1)*(1-zratio)
        write(8,"(2(',',g11.5),$)")pm1at,pm1ac

C****Starting Modeling Aerosol
         if(hourjulian.le.ichem1_hour_end) then
 	  aline2='CHEM3D'
	 else
	  aline2='CHEM3D2'
	 endif  
                  
	 do L=1,laero
 
         if(laexist(L).eq.1) then
 
	  if(.not.INTERP3(trim(aline2),name_aero(L),'pathway',
     1	   iyear*1000+jday,nowtime,imax*jmax*kchem3d,work)) then
           print*,'interp error for ',name_aero(L),' on ',jday,nowtime	 
	   stop
	  endif 
	 
!	 else          !    ichem1_hour_end  <  hourjulian< (ichem1_hour_end+ichemstep)
!	  idaytmp1=ichem1_hour_end/24
!	  itimetmp1=mod(ichem1_hour_end,24)
	  
!	  if(.not.READ3('CHEM3D',name_aero(L),ALLAYS3,
!     1     iyear*1000+idaytmp1,itimetmp1*10000,work)) stop      ! read last time step of the first file

!	  idaytmp2=(ichem1_hour_end+ichemstep)/24
!	  itimetmp2=mod(ichem1_hour_end+ichemstep,24)	 

!	  if(.not.READ3('CHEM3D2',name_aero(L),ALLAYS3,
!     1     iyear*1000+idaytmp2,itimetmp2*10000,work2)) stop      ! read the time step of the second file
	  
!	  tratio=(hourjulian-ichem1_hour_end)/ichemstep
!	  if(tratio.lt.0.or.tratio.gt.1) then
!	   print*,'wrong tratio',hourjulian,ichem1_hour_end,ichemstep
!	   stop
!	  endif 
!	  work(1:imax,1:jmax,1:kchem3d)=work(1:imax,1:jmax,1:kchem3d)
!     1      *(1-tratio) +  work2(1:imax,1:jmax,1:kmax)*tratio
 	 
	 do mlev=1,kchem3d
	  conlevs(mlev)=(1-yratio)*(work(int(x),int(y),mlev)*
     1	  (1-xratio)+work(int(x)+1,int(y),mlev)*xratio)+yratio*
     2    (work(int(x),int(y)+1,mlev)*(1-xratio)+
     3     work(int(x)+1,int(y)+1,mlev)*xratio)
         enddo

         conc(L)=conlevs(kp)*zratio+conlevs(kp-1)*(1-zratio)
         write(8,"(',',g11.5,$)")conc(L)*temp/press
     1	  *101300/273.0                               !convert ug/m3 to ug/std m3 

        else  ! laexist=0
	 conc(L)=0.
         write(8,"(',',$)")
        endif
         
         if(name_aero(L).eq.'AECJ') then
          BC=conc(L)*pm1ac+conc(L-1)*pm1at
          write(8,"(',',g11.5,$)")BC*temp/press
     1	  *101300/273.0                               !PM1 BC
         else if(name_aero(L).eq.'AOTHRJ') then
          AOTHR=conc(L)*pm1ac+conc(L-1)*pm1at
          write(8,"(',',g11.5,$)")AOTHR*temp/press
     1	  *101300/273.0                               !total Other PM2.5
	 else if(any((/'ACLK','ASO4K','ANO3K','ANAK',
     1	  'ANH4K','ASEACAT'/).eq.trim(name_aero(L)))) then
	  ftmp=conc(L-2)*pm1at+conc(L-1)*pm1ac        ! PM1 aerosol
	  if (name_aero(L).eq.'ASEACAT') then
	   ttmp=conc(L-1)+conc(L-2)+conc(L)*0.8373
	  else  
           ttmp=sum(conc(L-2:L))                       ! total aerosol
	  endif 
	  write(8,"(3(',',g11.5),$)")ftmp*temp/press
     1	  *101300/273.,ttmp*temp/press*101300/273.,
     2    ftmp/ttmp 
          if(name_aero(L).eq.'ANO3K') vnitrate=(ftmp+ctmp)*temp/    ! ug/m3 -> ppbv
     1	   press*101300/273.15*22.4/62.          ! NO3- molecular weight 
	 else if(any((/'ACAJ','AKJ','AFEJ','AALJ',
     1	  'ASIJ','ATIJ','AMNJ','AMGJ'/).eq.trim(name_aero(L)))) then
          write(8,"(',',g11.5,$)")pm1ac*conc(L)*temp/press*101300/273.
         endif   
	  
	 enddo   ! end of inorganic aerosol loop
         
         write(8,"(',',g11.5,$)")vnitrate+hno3    ! total nitrate (aerosol + gas)

!----POA and SOA species
         poai=0; soai=0; poaj=0; soaj=0	 
	 do L=1,lpomi+lsomi+lpomj+lsomj
 
	  if(.not.INTERP3(trim(aline2),name_tom(L),'pathway',
     1	   iyear*1000+jday,nowtime,imax*jmax*kchem3d,work)) stop
     	  do mlev=1,kchem3d
	   conlevs(mlev)=(1-yratio)*(work(int(x),int(y),mlev)*
     1	   (1-xratio)+work(int(x)+1,int(y),mlev)*xratio)+yratio*
     2     (work(int(x),int(y)+1,mlev)*(1-xratio)+
     3      work(int(x)+1,int(y)+1,mlev)*xratio)
          enddo
	  conc(L)=(conlevs(kp)*zratio+conlevs(kp-1)*(1-zratio))
	  if ( any(name_pomi.eq.name_tom(L))) then
	    poai=poai+conc(L)*temp/press*101300/273.0
	  else if ( any(name_somi.eq.name_tom(L))) then
	    soai=soai+conc(L)*temp/press*101300/273.0
	  else if ( any(name_pomj.eq.name_tom(L))) then
	    poaj=poaj+conc(L)*temp/press*101300/273.0
	  else if ( any(name_somj.eq.name_tom(L))) then
	    soaj=soaj+conc(L)*temp/press*101300/273.0
	  endif
	 enddo 
	  write(8,"(9(',',g11.5),$)")poai,soai,poai+soai,poaj,soaj,
     1	   poaj+soaj,poai*pm1at+poaj*pm1ac,soai*pm1at+soaj*pm1ac,
     2     (poai+soai)*pm1at+(poaj+soaj)*pm1ac
	    
	 
!! Jvalue
        if(usefile1) then
	  aline2='JVFILE'
	else
	  aline2='JVFILE2'
	endif
	do L=1,lj
	 if(.not.interp3(trim(aline2),name_j(L),'youhua',
     1	  iyear*1000+jday,nowtime,imax*jmax*kchem3d,work)) stop
     	  do mlev=1,kchem3d
	   conlevs(mlev)=(1-yratio)*(work(int(x),int(y),mlev)*
     1	   (1-xratio)+work(int(x)+1,int(y),mlev)*xratio)+yratio*
     2     (work(int(x),int(y)+1,mlev)*(1-xratio)+
     3      work(int(x)+1,int(y)+1,mlev)*xratio)
          enddo
	  conc(L)=(conlevs(kp)*zratio+conlevs(kp-1)*(1-zratio))
	  write(8,"(',',g11.5,$)")conc(L)/60.
	 enddo
!! AOP
        if(usefile1) then
          aline2='AOP'
        else
          aline2='AOP2'
        endif
        do L=1,laoe
         if(.not.interp3(trim(aline2),name_aoe(L),'youhua',
     1    iyear*1000+jday,nowtime,imax*jmax*kchem3d,work)) stop
          do mlev=1,kchem3d
           conlevs(mlev)=(1-yratio)*(work(int(x),int(y),mlev)*
     1     (1-xratio)+work(int(x)+1,int(y),mlev)*xratio)+yratio*
     2     (work(int(x),int(y)+1,mlev)*(1-xratio)+
     3      work(int(x)+1,int(y)+1,mlev)*xratio)
          enddo
          conc(L)=(conlevs(kp)*zratio+conlevs(kp-1)*(1-zratio))
          write(8,"(',',g11.5,$)")conc(L)
         enddo

 50	write(8,'(1x)')
	
       enddo                             ! end of one file 
 99    print*,'end of file ',afile(mn)
       close(7)
       close(8)

       print*,'end of total process'
       end




      REAL FUNCTION DAOD(MSS,rhu,lspecies)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  look-up-table for converting aerosol concentration to corresponding      =*
*=  optical extinction coffecient                                            =*
*-----------------------------------------------------------------------------*
*=  INPUT     :                                                              =*
*=  MSS     - REAL, aerosol concentration in ug/m3                           =*
*=  RHU     - REAL, Relative humidty (%)                                     =*
*=  lspecies - 1-dust, 2-sulfate, 3-black carbon, 4-sea salt                 =*
*=                        5-organice carbon                                  =*
*=  DAOD    - REAL, optical extinction coffecient                            =*   
*=            wavelength and aerosol type                                    =*
*-----------------------------------------------------------------------------*
*=  in order to calculate aerosol optical depth, require to input            =*
*=     aerosol mass concentration in (COL,ROW,each layer, aerosol type)      =*
*=     to MSS(aerosol type); and relative humidity for that grid cell.       =*
*-----------------------------------------------------------------------------*


*_________________________________________________
* altitude, wavelength grids size
      INTEGER KZ, KW
* altitude
      PARAMETER(KZ=101)
* wavelength
      PARAMETER(KW=130)

* aerosol layers, aerosol types
      INTEGER NAETP
* number of aerosol types 1=dust,2=watersoluble(50%sulfate),3=black carbon, 4=sea salt,
*                         5=organic carbon
      PARAMETER(NAETP=5)
      
      
*_________________________________________________

      
      REAL     MSS,rhu
      

* local:
      INTEGER  ierr
      INTEGER  IATP, ILAY, IRH, IW   !loop control

      REAL     RH00(8),  RHUM
      REAL     WAVEI(7)
      
      REAL     CONV(NAETP),NDS  ! aerosol numer density
      
      REAL     SBEA(6,8,NAETP), SBAA(6,8,NAETP), SGA(6,8,NAETP)
      REAL     ABEAR(7,NAETP), AGAR(7,NAETP), AOMR(7,NAETP)
      REAL     extaer(7), albaer(7), asmaer(7)


*_______________________________________________________________________
      DATA WAVEI /.185E+3,.25E+3,.3E+3,
     >               .4E+3,.55E+3,.7E+3,1.5E+3/
     
      DATA RH00  /.0,.50,.70,.80,.90,.95,.98,.99/

C---------------------------------------------------------------------------
C** CONV convert ug/m3 to #/cm3 for dust 1, sulfate 2, soot 3, sea salt 4 **
C** aerosol size distributions follow OPAC 
C** organic carbon 5, optical from Liousse et al., 1996 & Cooke et al., 1999
C** r=0.0212micron+-2.24, density=1.8g/cm3, optcial reliable between
C** 300~1060nm wavelength
C** dust    300.142 particles/cm3 / 221.8ug/m3 = 1.3532    (#/cm3)/(ug/m3)
C** sulfate 28000   particles/cm3 / 56.0 ug/m3 = 500.0     (#/cm3)/(ug/m3)
C** soot    130,000 particles/cm3 / 7.8  ug/m3 = 16666.667 (#/cm3)/(ug/m3)
C** seasalt 20.0032 particles/cm3 / 39.5 ug/m3 = 0.50641   (#/cm3)/(ug/m3)
C** organic carbon  no reliable source, use WaterSoluble data in OPAC
C**                 1.34E-3(micro-g/m3/part/cm3)=746.27    (#/cm3)/(ug/m3)

c      DATA CONV  / 1.3532E+00, 5.0000E+02, 1.6667E+04,
c     >             5.0641E-01, 7.4627E+02/
      DATA CONV  / 1.9000E+00, 5.0000E+02, 1.6667E+04,      ! increasing dust converting ratio
     >             5.0641E-01, 7.4627E+02/                  ! since we use it as fine dusts    

C---------------------------------------------------------------------------
C** SBEA extinction coefficient (1/km) normalized to 1 particle/cm3
      DATA SBEA  /
C** mineral dust                                 **
     >  4.153E-04,4.177E-04,4.190E-04,4.182E-04,4.180E-04,4.032E-04,
     >  4.153E-04,4.177E-04,4.190E-04,4.182E-04,4.180E-04,4.032E-04,
     >  4.153E-04,4.177E-04,4.190E-04,4.182E-04,4.180E-04,4.032E-04,
     >  4.153E-04,4.177E-04,4.190E-04,4.182E-04,4.180E-04,4.032E-04,
     >  4.153E-04,4.177E-04,4.190E-04,4.182E-04,4.180E-04,4.032E-04,
     >  4.153E-04,4.177E-04,4.190E-04,4.182E-04,4.180E-04,4.032E-04,
     >  4.153E-04,4.177E-04,4.190E-04,4.182E-04,4.180E-04,4.032E-04,
     >  4.153E-04,4.177E-04,4.190E-04,4.182E-04,4.180E-04,4.032E-04,
C** watersoluable particles  (50% sulfate)       **
     >  9.520E-06,8.182E-06,6.011E-06,3.905E-06,2.649E-06,5.937E-07,
     >  1.511E-05,1.301E-05,9.654E-06,6.366E-06,4.382E-06,9.935E-07,
     >  1.835E-05,1.587E-05,1.188E-05,7.913E-06,5.493E-06,1.264E-06,
     >  2.169E-05,1.885E-05,1.423E-05,9.584E-06,6.711E-06,1.571E-06,
     >  2.945E-05,2.589E-05,1.994E-05,1.374E-05,9.794E-06,2.392E-06,
     >  4.076E-05,3.637E-05,2.874E-05,2.038E-05,1.485E-05,3.846E-06,
     >  6.157E-05,5.616E-05,4.605E-05,3.404E-05,2.560E-05,7.265E-06,
     >  7.993E-05,7.410E-05,6.230E-05,4.740E-05,3.644E-05,1.105E-05,
C** black carbon                                 **
     >  2.2035E-06,1.97E-06,13.995E-07,9.0E-07,6.375E-07,2.491E-07,
     >  2.2035E-06,1.97E-06,13.995E-07,9.0E-07,6.375E-07,2.491E-07,
     >  2.2035E-06,1.97E-06,13.995E-07,9.0E-07,6.375E-07,2.491E-07,
     >  2.2035E-06,1.97E-06,13.995E-07,9.0E-07,6.375E-07,2.491E-07,
     >  2.2035E-06,1.97E-06,13.995E-07,9.0E-07,6.375E-07,2.491E-07,
     >  2.2035E-06,1.97E-06,13.995E-07,9.0E-07,6.375E-07,2.491E-07,
     >  2.2035E-06,1.97E-06,13.995E-07,9.0E-07,6.375E-07,2.491E-07,
     >  2.2035E-06,1.97E-06,13.995E-07,9.0E-07,6.375E-07,2.491E-07,
C** sea salt
     >  9.132E-04,9.455E-04,9.934E-04,1.037E-03,1.040E-03,7.527E-04,
     >  2.255E-03,2.304E-03,2.412E-03,2.536E-03,2.602E-03,2.213E-03,
     >  2.814E-03,2.880E-03,3.005E-03,3.164E-03,3.261E-03,2.919E-03,
     >  3.378E-03,3.445E-03,3.593E-03,3.777E-03,3.908E-03,3.643E-03,
     >  4.757E-03,4.825E-03,5.006E-03,5.252E-03,5.459E-03,5.441E-03,
     >  6.929E-03,7.044E-03,7.249E-03,7.560E-03,7.840E-03,8.297E-03,
     >  1.199E-02,1.214E-02,1.238E-02,1.280E-02,1.318E-02,1.462E-02,
     >  1.832E-02,1.849E-02,1.878E-02,1.930E-02,1.975E-02,2.214E-02,
C** organic carbon from Liousse et al., 1996, RH dependency follow
C** watersoluable particles.
     >  1.405E-05,11.72E-06,8.969E-06,6.700E-06,5.259E-06,2.541E-06,                    
     >  2.228E-05,1.861E-05,1.424E-05,10.63E-06,8.346E-06,4.033E-06,                    
     >  2.707E-05,2.260E-05,1.728E-05,1.291E-05,10.14E-06,4.897E-06,                    
     >  3.200E-05,2.138E-05,2.044E-05,1.526E-05,11.98E-06,5.790E-06,                    
     >  4.345E-05,3.627E-05,2.775E-05,2.072E-05,1.626E-05,7.861E-06,                    
     >  6.014E-05,5.020E-05,3.840E-05,2.869E-05,2.251E-05,10.88E-06,                    
     >  9.083E-05,7.582E-05,5.801E-05,4.333E-05,3.401E-05,1.643E-05,                    
     >  11.79E-05,9.844E-05,7.530E-05,5.625E-05,4.415E-05,2.133E-05/                    
     
      DATA SBAA  /
C** mineral dust                                 **
     >  6.234E-01,6.539E-01,7.549E-01,8.732E-01,9.125E-01,9.313E-01,
     >  6.234E-01,6.539E-01,7.549E-01,8.732E-01,9.125E-01,9.313E-01,
     >  6.234E-01,6.539E-01,7.549E-01,8.732E-01,9.125E-01,9.313E-01,
     >  6.234E-01,6.539E-01,7.549E-01,8.732E-01,9.125E-01,9.313E-01,
     >  6.234E-01,6.539E-01,7.549E-01,8.732E-01,9.125E-01,9.313E-01,
     >  6.234E-01,6.539E-01,7.549E-01,8.732E-01,9.125E-01,9.313E-01,
     >  6.234E-01,6.539E-01,7.549E-01,8.732E-01,9.125E-01,9.313E-01,
     >  6.234E-01,6.539E-01,7.549E-01,8.732E-01,9.125E-01,9.313E-01,
C** watersoluable particles  (50% sulfate)       **
     >  8.360E-01,9.493E-01,9.685E-01,9.615E-01,9.522E-01,7.631E-01,
     >  8.918E-01,9.687E-01,9.808E-01,9.765E-01,9.708E-01,8.479E-01,
     >  9.094E-01,9.744E-01,9.844E-01,9.811E-01,9.765E-01,8.773E-01,
     >  9.223E-01,9.785E-01,9.871E-01,9.843E-01,9.806E-01,8.990E-01,
     >  9.415E-01,9.842E-01,9.907E-01,9.890E-01,9.865E-01,9.308E-01,
     >  9.568E-01,9.888E-01,9.935E-01,9.925E-01,9.910E-01,9.548E-01,
     >  9.707E-01,9.926E-01,9.959E-01,9.954E-01,9.946E-01,9.741E-01,
     >  9.773E-01,9.944E-01,9.970E-01,9.967E-01,9.962E-01,9.819E-01,
C** black carbon                                 **
     >  3.081E-01,3.128E-01,2.672E-01,2.088E-01,1.616E-01,4.407E-02,
     >  3.081E-01,3.128E-01,2.672E-01,2.088E-01,1.616E-01,4.407E-02,
     >  3.081E-01,3.128E-01,2.672E-01,2.088E-01,1.616E-01,4.407E-02,
     >  3.081E-01,3.128E-01,2.672E-01,2.088E-01,1.616E-01,4.407E-02,
     >  3.081E-01,3.128E-01,2.672E-01,2.088E-01,1.616E-01,4.407E-02,
     >  3.081E-01,3.128E-01,2.672E-01,2.088E-01,1.616E-01,4.407E-02,
     >  3.081E-01,3.128E-01,2.672E-01,2.088E-01,1.616E-01,4.407E-02,
     >  3.081E-01,3.128E-01,2.672E-01,2.088E-01,1.616E-01,4.407E-02,
C ** sea salt
     >  9.998E-01,9.999E-01,1.000E+00,1.000E+00,1.000E+00,9.956E-01,
     >  9.999E-01,1.000E+00,1.000E+00,1.000E+00,1.000E+00,9.970E-01,
     >  9.999E-01,1.000E+00,1.000E+00,1.000E+00,1.000E+00,9.971E-01,
     >  1.000E+00,1.000E+00,1.000E+00,1.000E+00,1.000E+00,9.971E-01,
     >  1.000E+00,1.000E+00,1.000E+00,1.000E+00,1.000E+00,9.969E-01,
     >  1.000E+00,1.000E+00,1.000E+00,1.000E+00,1.000E+00,9.966E-01,
     >  1.000E+00,1.000E+00,1.000E+00,1.000E+00,1.000E+00,9.957E-01,
     >  1.000E+00,1.000E+00,1.000E+00,1.000E+00,1.000E+00,9.947E-01,
C ** organic carbon from Cooke et al., 1999 only 550nm is available,
C ** at other wavelength and RH, is the same as 550nm, be careful.
     >  9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,
     >  9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,
     >  9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,
     >  9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,
     >  9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,
     >  9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,
     >  9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,
     >  9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01/
     
     
      DATA SGA   /
C** mineral dust                                 **
     >  8.572E-01,8.363E-01,7.849E-01,7.329E-01,7.112E-01,6.933E-01,
     >  8.572E-01,8.363E-01,7.849E-01,7.329E-01,7.112E-01,6.933E-01,
     >  8.572E-01,8.363E-01,7.849E-01,7.329E-01,7.112E-01,6.933E-01,
     >  8.572E-01,8.363E-01,7.849E-01,7.329E-01,7.112E-01,6.933E-01,
     >  8.572E-01,8.363E-01,7.849E-01,7.329E-01,7.112E-01,6.933E-01,
     >  8.572E-01,8.363E-01,7.849E-01,7.329E-01,7.112E-01,6.933E-01,
     >  8.572E-01,8.363E-01,7.849E-01,7.329E-01,7.112E-01,6.933E-01,
     >  8.572E-01,8.363E-01,7.849E-01,7.329E-01,7.112E-01,6.933E-01,
C** watersoluable particles  (50% sulfate)       **
     >  6.900E-01,6.590E-01,6.390E-01,6.140E-01,5.890E-01,4.850E-01,
     >  7.290E-01,7.110E-01,6.950E-01,6.720E-01,6.470E-01,5.380E-01,
     >  7.400E-01,7.260E-01,7.120E-01,6.900E-01,6.660E-01,5.590E-01,
     >  7.480E-01,7.360E-01,7.240E-01,7.040E-01,6.810E-01,5.750E-01,
     >  7.590E-01,7.520E-01,7.430E-01,7.250E-01,7.040E-01,6.030E-01,
     >  7.670E-01,7.640E-01,7.580E-01,7.430E-01,7.240E-01,6.310E-01,
     >  7.730E-01,7.740E-01,7.720E-01,7.610E-01,7.450E-01,6.630E-01,
     >  7.760E-01,7.780E-01,7.780E-01,7.700E-01,7.560E-01,6.820E-01,
C** black carbon                                 **
     >  5.020E-01,4.530E-01,3.960E-01,3.360E-01,2.910E-01,1.580E-01,
     >  5.020E-01,4.530E-01,3.960E-01,3.360E-01,2.910E-01,1.580E-01,
     >  5.020E-01,4.530E-01,3.960E-01,3.360E-01,2.910E-01,1.580E-01,
     >  5.020E-01,4.530E-01,3.960E-01,3.360E-01,2.910E-01,1.580E-01,
     >  5.020E-01,4.530E-01,3.960E-01,3.360E-01,2.910E-01,1.580E-01,
     >  5.020E-01,4.530E-01,3.960E-01,3.360E-01,2.910E-01,1.580E-01,
     >  5.020E-01,4.530E-01,3.960E-01,3.360E-01,2.910E-01,1.580E-01,
     >  5.020E-01,4.530E-01,3.960E-01,3.360E-01,2.910E-01,1.580E-01,
C** sea salt
     >  7.219E-01,7.070E-01,6.990E-01,6.919E-01,6.968E-01,7.058E-01,
     >  7.856E-01,7.826E-01,7.757E-01,7.717E-01,7.736E-01,7.844E-01,
     >  7.985E-01,7.906E-01,7.837E-01,7.787E-01,7.806E-01,7.924E-01,
     >  8.034E-01,7.996E-01,7.896E-01,7.847E-01,7.856E-01,7.974E-01,
     >  8.163E-01,8.124E-01,8.006E-01,7.936E-01,7.906E-01,8.025E-01,
     >  8.272E-01,8.223E-01,8.115E-01,8.016E-01,7.976E-01,8.055E-01,
     >  8.380E-01,8.361E-01,8.283E-01,8.146E-01,8.066E-01,8.056E-01,
     >  8.447E-01,8.439E-01,8.382E-01,8.255E-01,8.156E-01,8.056E-01,
C** organic carbon no assymmetry factor is unavailable, use g value
C** for continental clean at 550nm from OPAC, then independent of
C** wavlength and RH, can not be used for radiative transfer calculation.
     >  7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,
     >  7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,
     >  7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,
     >  7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,
     >  7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,
     >  7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,
     >  7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,
     >  7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01/
      
C-------------------------------------------------------------------------C
C   SBEA(6,8,5) extinction coefficient                                    C
C   SBAA(6,8,5) single scattering albedo                                  C
C   SGA (6,8,5) asymmetry factor                                          C
C        6 wavelength, .25, .3, .4, .55, .7, 1.5 micron                   C
C        8 relative humidity, .0, .5, .7, .8, .9, .95, .98, .99           C
C        5 aerosol type, 1-dust, 2-sulfate, 3-black carbon, 4-sea salt    C
C                        5-organice carbon                                C
C-------------------------------------------------------------------------C
C   optical parameters are calculated using OPAC                          **
C ** 1) mineral 300.142 particles/cm3, 221.8 ug/m3                        **
C **    number density: nuc 269.5, acc 30.5, coa 0.142 /cm3               **
C **    mass %: nuc 7.5 (3.4%), acc 168.7 (76.1%), coa 45.6 (20.5%) ug/m3 **
C ** 2) water soluble  28000 particles/cm3, 56.0 ug/m3                    **
C ** 3) black carbon 130,000 particles/cm3,  7.8 ug/m3                    **
C ** 4) sea salt 20.0032 part/cm3, 39.5ug/m3                              **
C **    number density: acc 20.(99.98%), coa 3.2E-3(0.02%)                **
C **    mass : acc 38.6 (97.72%), coa 0.9 (2.28%)                         ** 
C ** 5) organic carbon as Watersoluable in OPAC,1.34ug/m3/part/cm3        **
C-------------------------------------------------------------------------C



C--- aerosol mass concentration input here in unit of micro-g/m3 ---------C
C    IATP=1, dust;          IATP=2, watersoluble (50% sulfate);           C
C    IATP=3, black carbon;  IATP=4, sea salt;                             C
C    For example:                                                         C
C    MSS(ILAY,1) = total dust concentration at layer ILAY                 C
C--- aerosol mass concentration input END --------------------------------C
        
C--- convert aerosol mass ug/m3 to #/cm3 ---------------------------------C
        NDS =MSS *CONV(lspecies)
     

C--- get relative humidity value -----------------------------------------C
C    RH value rhu(COL,ROW,ILAY) is required                               C
C      rhu(ILAY)= 0.
      
      RHUM = rhu/100.
      
      IF(RHUM.GT.0.99) RHUM=0.99
      IF(RHUM.LT.0.00) RHUM=0.00
      
      IATP=lspecies
       
       DO IW =2, 7                   !wavelength loop
         ABEAR(IW,IATP)=0.
         AOMR (IW,IATP)=0.
         AGAR (IW,IATP)=0.

       IF(IATP.EQ.2.OR.IATP.EQ.4) THEN        !RH sensitive
        DO IRH=1,7                            !relative humidity loop
        IF((RHUM.GE.RH00(IRH)).AND.(RHUM.LE.RH00(IRH+1))) THEN
         ABEAR(IW,IATP)=SBEA(IW-1,IRH,IATP)+ ((RHUM-RH00(IRH))/
     >            (RH00(IRH+1)-RH00(IRH)))*
     >            (SBEA(IW-1,IRH+1,IATP)-SBEA(IW-1,IRH,IATP))
         AOMR (IW,IATP)=SBAA(IW-1,IRH,IATP)+ ((RHUM-RH00(IRH))/
     >            (RH00(IRH+1)-RH00(IRH)))*
     >            (SBAA(IW-1,IRH+1,IATP)-SBAA(IW-1,IRH,IATP))
         AGAR (IW,IATP)=SGA (IW-1,IRH,IATP)+ ((RHUM-RH00(IRH))/
     >            (RH00(IRH+1)-RH00(IRH)))*
     >            (SGA (IW-1,IRH+1,IATP)-SGA (IW-1,IRH,IATP))
         END IF
        END DO                                !relative humidity loop END
       ELSE                                   !RH independent
         ABEAR(IW,IATP)=SBEA(IW-1,1,IATP)
         AOMR (IW,IATP)=SBAA(IW-1,1,IATP)
         AGAR (IW,IATP)=SGA (IW-1,1,IATP)
       END IF                                 !RH ENDIF
        
       END DO                        !wavelength loop END

       ABEAR(1,IATP)= ABEAR(2,IATP)+(ABEAR(2,IATP)-ABEAR(3,IATP))*
     >          (wavei(1)-wavei(2))/(wavei(2)-wavei(3))
       AOMR (1,IATP)= AOMR (2,IATP)+(AOMR (2,IATP) -AOMR(3,IATP))*
     >          (wavei(1)-wavei(2))/(wavei(2)-wavei(3))
       AGAR (1,IATP)= AGAR (2,IATP)+(AGAR (2,IATP) -AGAR(3,IATP))*
     >          (wavei(1)-wavei(2))/(wavei(2)-wavei(3))
     
       DO IW=1,7                     !wavelength loop
        extaer(IW)=ABEAR(IW,lspecies)
        albaer(IW)=AOMR (IW,lspecies)
        asmaer(IW)=AGAR (IW,lspecies)        
       END DO                        !wavelength loop END
       
      do IW=2,7
C---- DAOD is at 550nm wavelength       
      if(wavei(iw-1).le.550.and.wavei(iw).ge.550) then
       daod=(extaer(iw-1)+(extaer(iw)-extaer(iw-1))/
     1  (wavei(iw)-wavei(iw-1))*(550-wavei(iw-1)))*NDS
       return 
      endif
      enddo 
      END

      subroutine obsscat(MSS,rhu,tabsorp,tscatter)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  look-up-table for converting aerosol concentration to total              =*
*=  single scatter albedo                                                    =*
*-----------------------------------------------------------------------------*
*=  INPUT     :                                                              =*
*=  MSS     - REAL, aerosol concentration in ug/m3                           =*
*=  RHU     - REAL, Relative humidty (%)                                     =*
*=  lspecies - 1-dust, 2-sulfate, 3-black carbon, 4-sea salt                 =*
*=                        5-organice carbon                                  =*
*=  salbedo    - REAL, optical albedo coefficient                            =*
*-----------------------------------------------------------------------------*
*=  in order to calculate aerosol optical depth, require to input            =*
*=     aerosol mass concentration in (COL,ROW,each layer, aerosol type)      =*
*=     to MSS(aerosol type); and relative humidity for that grid cell.       =*
*-----------------------------------------------------------------------------*


*_________________________________________________
* altitude, wavelength grids size
      INTEGER KZ, KW
* altitude
      PARAMETER(KZ=101)
* wavelength
      PARAMETER(KW=130)

* aerosol layers, aerosol types
      INTEGER NAETP
* number of aerosol types 1=dust,2=watersoluble(50%sulfate),3=black carbon, 4=sea salt,
*                         5=organic carbon
      PARAMETER(NAETP=5)
      
      
*_________________________________________________

      
      REAL     MSS(NAETP),rhu
      

* local:
      INTEGER  ierr
      INTEGER  IATP, ILAY, IRH, IW   !loop control

      REAL     RH00(8),  RHUM
      REAL     WAVEI(7)
      
      REAL     CONV(NAETP),NDS(NAETP), NAER  ! aerosol numer density
      
      REAL     SBEA(6,8,NAETP), SBAA(6,8,NAETP), SGA(6,8,NAETP)
      REAL     ABEAR(7,NAETP), AGAR(7,NAETP), AOMR(7,NAETP)
      REAL     extaer(7), albaer(7), asmaer(7), absorp(7)


*_______________________________________________________________________
      DATA WAVEI /.185E+3,.25E+3,.3E+3,
     >               .4E+3,.55E+3,.7E+3,1.5E+3/
     
      DATA RH00  /.0,.50,.70,.80,.90,.95,.98,.99/

C---------------------------------------------------------------------------
C** CONV convert ug/m3 to #/cm3 for dust 1, sulfate 2, soot 3, sea salt 4 **
C** aerosol size distributions follow OPAC 
C** organic carbon 5, optical from Liousse et al., 1996 & Cooke et al., 1999
C** r=0.0212micron+-2.24, density=1.8g/cm3, optcial reliable between
C** 300~1060nm wavelength
C** dust    300.142 particles/cm3 / 221.8ug/m3 = 1.3532    (#/cm3)/(ug/m3)
C** sulfate 28000   particles/cm3 / 56.0 ug/m3 = 500.0     (#/cm3)/(ug/m3)
C** soot    130,000 particles/cm3 / 7.8  ug/m3 = 16666.667 (#/cm3)/(ug/m3)
C** seasalt 20.0032 particles/cm3 / 39.5 ug/m3 = 0.50641   (#/cm3)/(ug/m3)
C** organic carbon  no reliable source, use WaterSoluble data in OPAC
C**                 1.34E-3(micro-g/m3/part/cm3)=746.27    (#/cm3)/(ug/m3)

c      DATA CONV  / 1.3532E+00, 5.0000E+02, 1.6667E+04,
c     >             5.0641E-01, 7.4627E+02/
      DATA CONV  / 1.9000E+00, 5.0000E+02, 1.6667E+04,      ! increasing dust converting ratio
     >             5.0641E-01, 7.4627E+02/                  ! since we use it as fine dusts    

C---------------------------------------------------------------------------
C** SBEA extinction coefficient (1/km) normalized to 1 particle/cm3
      DATA SBEA  /
C** mineral dust                                 **
     >  4.153E-04,4.177E-04,4.190E-04,4.182E-04,4.180E-04,4.032E-04,
     >  4.153E-04,4.177E-04,4.190E-04,4.182E-04,4.180E-04,4.032E-04,
     >  4.153E-04,4.177E-04,4.190E-04,4.182E-04,4.180E-04,4.032E-04,
     >  4.153E-04,4.177E-04,4.190E-04,4.182E-04,4.180E-04,4.032E-04,
     >  4.153E-04,4.177E-04,4.190E-04,4.182E-04,4.180E-04,4.032E-04,
     >  4.153E-04,4.177E-04,4.190E-04,4.182E-04,4.180E-04,4.032E-04,
     >  4.153E-04,4.177E-04,4.190E-04,4.182E-04,4.180E-04,4.032E-04,
     >  4.153E-04,4.177E-04,4.190E-04,4.182E-04,4.180E-04,4.032E-04,
C** watersoluable particles  (50% sulfate)       **
     >  9.520E-06,8.182E-06,6.011E-06,3.905E-06,2.649E-06,5.937E-07,
     >  1.511E-05,1.301E-05,9.654E-06,6.366E-06,4.382E-06,9.935E-07,
     >  1.835E-05,1.587E-05,1.188E-05,7.913E-06,5.493E-06,1.264E-06,
     >  2.169E-05,1.885E-05,1.423E-05,9.584E-06,6.711E-06,1.571E-06,
     >  2.945E-05,2.589E-05,1.994E-05,1.374E-05,9.794E-06,2.392E-06,
     >  4.076E-05,3.637E-05,2.874E-05,2.038E-05,1.485E-05,3.846E-06,
     >  6.157E-05,5.616E-05,4.605E-05,3.404E-05,2.560E-05,7.265E-06,
     >  7.993E-05,7.410E-05,6.230E-05,4.740E-05,3.644E-05,1.105E-05,
C** black carbon                                 **
     >  2.2035E-06,1.97E-06,13.995E-07,9.0E-07,6.375E-07,2.491E-07,
     >  2.2035E-06,1.97E-06,13.995E-07,9.0E-07,6.375E-07,2.491E-07,
     >  2.2035E-06,1.97E-06,13.995E-07,9.0E-07,6.375E-07,2.491E-07,
     >  2.2035E-06,1.97E-06,13.995E-07,9.0E-07,6.375E-07,2.491E-07,
     >  2.2035E-06,1.97E-06,13.995E-07,9.0E-07,6.375E-07,2.491E-07,
     >  2.2035E-06,1.97E-06,13.995E-07,9.0E-07,6.375E-07,2.491E-07,
     >  2.2035E-06,1.97E-06,13.995E-07,9.0E-07,6.375E-07,2.491E-07,
     >  2.2035E-06,1.97E-06,13.995E-07,9.0E-07,6.375E-07,2.491E-07,
C** sea salt
     >  9.132E-04,9.455E-04,9.934E-04,1.037E-03,1.040E-03,7.527E-04,
     >  2.255E-03,2.304E-03,2.412E-03,2.536E-03,2.602E-03,2.213E-03,
     >  2.814E-03,2.880E-03,3.005E-03,3.164E-03,3.261E-03,2.919E-03,
     >  3.378E-03,3.445E-03,3.593E-03,3.777E-03,3.908E-03,3.643E-03,
     >  4.757E-03,4.825E-03,5.006E-03,5.252E-03,5.459E-03,5.441E-03,
     >  6.929E-03,7.044E-03,7.249E-03,7.560E-03,7.840E-03,8.297E-03,
     >  1.199E-02,1.214E-02,1.238E-02,1.280E-02,1.318E-02,1.462E-02,
     >  1.832E-02,1.849E-02,1.878E-02,1.930E-02,1.975E-02,2.214E-02,
C** organic carbon from Liousse et al., 1996, RH dependency follow
C** watersoluable particles.
     >  1.405E-05,11.72E-06,8.969E-06,6.700E-06,5.259E-06,2.541E-06,                    
     >  2.228E-05,1.861E-05,1.424E-05,10.63E-06,8.346E-06,4.033E-06,                    
     >  2.707E-05,2.260E-05,1.728E-05,1.291E-05,10.14E-06,4.897E-06,                    
     >  3.200E-05,2.138E-05,2.044E-05,1.526E-05,11.98E-06,5.790E-06,                    
     >  4.345E-05,3.627E-05,2.775E-05,2.072E-05,1.626E-05,7.861E-06,                    
     >  6.014E-05,5.020E-05,3.840E-05,2.869E-05,2.251E-05,10.88E-06,                    
     >  9.083E-05,7.582E-05,5.801E-05,4.333E-05,3.401E-05,1.643E-05,                    
     >  11.79E-05,9.844E-05,7.530E-05,5.625E-05,4.415E-05,2.133E-05/                    
     
      DATA SBAA  /
C** mineral dust                                 **
     >  6.234E-01,6.539E-01,7.549E-01,8.732E-01,9.125E-01,9.313E-01,
     >  6.234E-01,6.539E-01,7.549E-01,8.732E-01,9.125E-01,9.313E-01,
     >  6.234E-01,6.539E-01,7.549E-01,8.732E-01,9.125E-01,9.313E-01,
     >  6.234E-01,6.539E-01,7.549E-01,8.732E-01,9.125E-01,9.313E-01,
     >  6.234E-01,6.539E-01,7.549E-01,8.732E-01,9.125E-01,9.313E-01,
     >  6.234E-01,6.539E-01,7.549E-01,8.732E-01,9.125E-01,9.313E-01,
     >  6.234E-01,6.539E-01,7.549E-01,8.732E-01,9.125E-01,9.313E-01,
     >  6.234E-01,6.539E-01,7.549E-01,8.732E-01,9.125E-01,9.313E-01,
C** watersoluable particles  (50% sulfate)       **
     >  8.360E-01,9.493E-01,9.685E-01,9.615E-01,9.522E-01,7.631E-01,
     >  8.918E-01,9.687E-01,9.808E-01,9.765E-01,9.708E-01,8.479E-01,
     >  9.094E-01,9.744E-01,9.844E-01,9.811E-01,9.765E-01,8.773E-01,
     >  9.223E-01,9.785E-01,9.871E-01,9.843E-01,9.806E-01,8.990E-01,
     >  9.415E-01,9.842E-01,9.907E-01,9.890E-01,9.865E-01,9.308E-01,
     >  9.568E-01,9.888E-01,9.935E-01,9.925E-01,9.910E-01,9.548E-01,
     >  9.707E-01,9.926E-01,9.959E-01,9.954E-01,9.946E-01,9.741E-01,
     >  9.773E-01,9.944E-01,9.970E-01,9.967E-01,9.962E-01,9.819E-01,
C** black carbon                                 **
     >  3.081E-01,3.128E-01,2.672E-01,2.088E-01,1.616E-01,4.407E-02,
     >  3.081E-01,3.128E-01,2.672E-01,2.088E-01,1.616E-01,4.407E-02,
     >  3.081E-01,3.128E-01,2.672E-01,2.088E-01,1.616E-01,4.407E-02,
     >  3.081E-01,3.128E-01,2.672E-01,2.088E-01,1.616E-01,4.407E-02,
     >  3.081E-01,3.128E-01,2.672E-01,2.088E-01,1.616E-01,4.407E-02,
     >  3.081E-01,3.128E-01,2.672E-01,2.088E-01,1.616E-01,4.407E-02,
     >  3.081E-01,3.128E-01,2.672E-01,2.088E-01,1.616E-01,4.407E-02,
     >  3.081E-01,3.128E-01,2.672E-01,2.088E-01,1.616E-01,4.407E-02,
C ** sea salt
     >  9.998E-01,9.999E-01,1.000E+00,1.000E+00,1.000E+00,9.956E-01,
     >  9.999E-01,1.000E+00,1.000E+00,1.000E+00,1.000E+00,9.970E-01,
     >  9.999E-01,1.000E+00,1.000E+00,1.000E+00,1.000E+00,9.971E-01,
     >  1.000E+00,1.000E+00,1.000E+00,1.000E+00,1.000E+00,9.971E-01,
     >  1.000E+00,1.000E+00,1.000E+00,1.000E+00,1.000E+00,9.969E-01,
     >  1.000E+00,1.000E+00,1.000E+00,1.000E+00,1.000E+00,9.966E-01,
     >  1.000E+00,1.000E+00,1.000E+00,1.000E+00,1.000E+00,9.957E-01,
     >  1.000E+00,1.000E+00,1.000E+00,1.000E+00,1.000E+00,9.947E-01,
C ** organic carbon from Cooke et al., 1999 only 550nm is available,
C ** at other wavelength and RH, is the same as 550nm, be careful.
     >  9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,
     >  9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,
     >  9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,
     >  9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,
     >  9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,
     >  9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,
     >  9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,
     >  9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01/
     
     
      DATA SGA   /
C** mineral dust                                 **
     >  8.572E-01,8.363E-01,7.849E-01,7.329E-01,7.112E-01,6.933E-01,
     >  8.572E-01,8.363E-01,7.849E-01,7.329E-01,7.112E-01,6.933E-01,
     >  8.572E-01,8.363E-01,7.849E-01,7.329E-01,7.112E-01,6.933E-01,
     >  8.572E-01,8.363E-01,7.849E-01,7.329E-01,7.112E-01,6.933E-01,
     >  8.572E-01,8.363E-01,7.849E-01,7.329E-01,7.112E-01,6.933E-01,
     >  8.572E-01,8.363E-01,7.849E-01,7.329E-01,7.112E-01,6.933E-01,
     >  8.572E-01,8.363E-01,7.849E-01,7.329E-01,7.112E-01,6.933E-01,
     >  8.572E-01,8.363E-01,7.849E-01,7.329E-01,7.112E-01,6.933E-01,
C** watersoluable particles  (50% sulfate)       **
     >  6.900E-01,6.590E-01,6.390E-01,6.140E-01,5.890E-01,4.850E-01,
     >  7.290E-01,7.110E-01,6.950E-01,6.720E-01,6.470E-01,5.380E-01,
     >  7.400E-01,7.260E-01,7.120E-01,6.900E-01,6.660E-01,5.590E-01,
     >  7.480E-01,7.360E-01,7.240E-01,7.040E-01,6.810E-01,5.750E-01,
     >  7.590E-01,7.520E-01,7.430E-01,7.250E-01,7.040E-01,6.030E-01,
     >  7.670E-01,7.640E-01,7.580E-01,7.430E-01,7.240E-01,6.310E-01,
     >  7.730E-01,7.740E-01,7.720E-01,7.610E-01,7.450E-01,6.630E-01,
     >  7.760E-01,7.780E-01,7.780E-01,7.700E-01,7.560E-01,6.820E-01,
C** black carbon                                 **
     >  5.020E-01,4.530E-01,3.960E-01,3.360E-01,2.910E-01,1.580E-01,
     >  5.020E-01,4.530E-01,3.960E-01,3.360E-01,2.910E-01,1.580E-01,
     >  5.020E-01,4.530E-01,3.960E-01,3.360E-01,2.910E-01,1.580E-01,
     >  5.020E-01,4.530E-01,3.960E-01,3.360E-01,2.910E-01,1.580E-01,
     >  5.020E-01,4.530E-01,3.960E-01,3.360E-01,2.910E-01,1.580E-01,
     >  5.020E-01,4.530E-01,3.960E-01,3.360E-01,2.910E-01,1.580E-01,
     >  5.020E-01,4.530E-01,3.960E-01,3.360E-01,2.910E-01,1.580E-01,
     >  5.020E-01,4.530E-01,3.960E-01,3.360E-01,2.910E-01,1.580E-01,
C** sea salt
     >  7.219E-01,7.070E-01,6.990E-01,6.919E-01,6.968E-01,7.058E-01,
     >  7.856E-01,7.826E-01,7.757E-01,7.717E-01,7.736E-01,7.844E-01,
     >  7.985E-01,7.906E-01,7.837E-01,7.787E-01,7.806E-01,7.924E-01,
     >  8.034E-01,7.996E-01,7.896E-01,7.847E-01,7.856E-01,7.974E-01,
     >  8.163E-01,8.124E-01,8.006E-01,7.936E-01,7.906E-01,8.025E-01,
     >  8.272E-01,8.223E-01,8.115E-01,8.016E-01,7.976E-01,8.055E-01,
     >  8.380E-01,8.361E-01,8.283E-01,8.146E-01,8.066E-01,8.056E-01,
     >  8.447E-01,8.439E-01,8.382E-01,8.255E-01,8.156E-01,8.056E-01,
C** organic carbon no assymmetry factor is unavailable, use g value
C** for continental clean at 550nm from OPAC, then independent of
C** wavlength and RH, can not be used for radiative transfer calculation.
     >  7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,
     >  7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,
     >  7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,
     >  7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,
     >  7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,
     >  7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,
     >  7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,
     >  7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01/
      
C-------------------------------------------------------------------------C
C   SBEA(6,8,5) extinction coefficient                                    C
C   SBAA(6,8,5) single scattering albedo                                  C
C   SGA (6,8,5) asymmetry factor                                          C
C        6 wavelength, .25, .3, .4, .55, .7, 1.5 micron                   C
C        8 relative humidity, .0, .5, .7, .8, .9, .95, .98, .99           C
C        5 aerosol type, 1-dust, 2-sulfate, 3-black carbon, 4-sea salt    C
C                        5-organice carbon                                C
C-------------------------------------------------------------------------C
C   optical parameters are calculated using OPAC                          **
C ** 1) mineral 300.142 particles/cm3, 221.8 ug/m3                        **
C **    number density: nuc 269.5, acc 30.5, coa 0.142 /cm3               **
C **    mass %: nuc 7.5 (3.4%), acc 168.7 (76.1%), coa 45.6 (20.5%) ug/m3 **
C ** 2) water soluble  28000 particles/cm3, 56.0 ug/m3                    **
C ** 3) black carbon 130,000 particles/cm3,  7.8 ug/m3                    **
C ** 4) sea salt 20.0032 part/cm3, 39.5ug/m3                              **
C **    number density: acc 20.(99.98%), coa 3.2E-3(0.02%)                **
C **    mass : acc 38.6 (97.72%), coa 0.9 (2.28%)                         ** 
C ** 5) organic carbon as Watersoluable in OPAC,1.34ug/m3/part/cm3        **
C-------------------------------------------------------------------------C



C--- aerosol mass concentration input here in unit of micro-g/m3 ---------C
C    IATP=1, dust;          IATP=2, watersoluble (50% sulfate);           C
C    IATP=3, black carbon;  IATP=4, sea salt;                             C
C    For example:                                                         C
C    MSS(ILAY,1) = total dust concentration at layer ILAY                 C
C--- aerosol mass concentration input END --------------------------------C
        
C--- convert aerosol mass ug/m3 to #/cm3 ---------------------------------C
      NDS(1:NAETP) =MSS(1:NAETP) *CONV(1:NAETP)
      
      Naer=0.
      DO IATP=1,NAETP
       NAER=NAER+NDS(IATP)
      ENDDO
     
C--- get relative humidity value -----------------------------------------C
C    RH value rhu(COL,ROW,ILAY) is required                               C
C      rhu(ILAY)= 0.
      
      RHUM = rhu/100.
      
      IF(RHUM.GT.0.99) RHUM=0.99
      IF(RHUM.LT.0.00) RHUM=0.00
      

      DO IATP=1, NAETP               !aerosol type loop
         ABEAR(1,IATP)=0.
         AOMR (1,IATP)=0.
         AGAR (1,IATP)=0.
      
       
       DO IW =2, 7                   !wavelength loop
         ABEAR(IW,IATP)=0.
         AOMR (IW,IATP)=0.
         AGAR (IW,IATP)=0.

       IF(IATP.EQ.2.OR.IATP.EQ.4) THEN        !RH sensitive
        DO IRH=1,7                            !relative humidity loop
        IF((RHUM.GE.RH00(IRH)).AND.(RHUM.LE.RH00(IRH+1))) THEN
         ABEAR(IW,IATP)=SBEA(IW-1,IRH,IATP)+ ((RHUM-RH00(IRH))/
     >            (RH00(IRH+1)-RH00(IRH)))*
     >            (SBEA(IW-1,IRH+1,IATP)-SBEA(IW-1,IRH,IATP))
         AOMR (IW,IATP)=SBAA(IW-1,IRH,IATP)+ ((RHUM-RH00(IRH))/
     >            (RH00(IRH+1)-RH00(IRH)))*
     >            (SBAA(IW-1,IRH+1,IATP)-SBAA(IW-1,IRH,IATP))
         AGAR (IW,IATP)=SGA (IW-1,IRH,IATP)+ ((RHUM-RH00(IRH))/
     >            (RH00(IRH+1)-RH00(IRH)))*
     >            (SGA (IW-1,IRH+1,IATP)-SGA (IW-1,IRH,IATP))
         END IF
        END DO                                !relative humidity loop END
       ELSE                                   !RH independent
         ABEAR(IW,IATP)=SBEA(IW-1,1,IATP)
         AOMR (IW,IATP)=SBAA(IW-1,1,IATP)
         AGAR (IW,IATP)=SGA (IW-1,1,IATP)
       END IF                                 !RH ENDIF
        
       END DO                        !wavelength loop END

       ABEAR(1,IATP)= ABEAR(2,IATP)+(ABEAR(2,IATP)-ABEAR(3,IATP))*
     >          (wavei(1)-wavei(2))/(wavei(2)-wavei(3))
       AOMR (1,IATP)= AOMR (2,IATP)+(AOMR (2,IATP) -AOMR(3,IATP))*
     >          (wavei(1)-wavei(2))/(wavei(2)-wavei(3))
       AGAR (1,IATP)= AGAR (2,IATP)+(AGAR (2,IATP) -AGAR(3,IATP))*
     >          (wavei(1)-wavei(2))/(wavei(2)-wavei(3))
       END DO                         !aerosol type loop END
     
      IF(Naer.GE.1.E-6) THEN        !IF Naer, avoid zero
C--- aerosol external mixture --------------------------------------------C
       DO IW=1,7                    !wavelength loop
        extaer(IW)=0.
	albaer(IW)=0.
	asmaer(IW)=0.
	absorp(IW)=0.
	
        DO IATP=1, NAETP                                      !aerosol type loop
         extaer(IW)= extaer(IW) + ABEAR(IW,IATP)*NDS(IATP)    ! total aerosol extinction  
         albaer(IW)= albaer(IW) + AOMR(IW,IATP)*NDS(IATP)
	 absorp(IW)= absorp(IW) + ABEAR(IW,IATP)*NDS(IATP)*   ! absorption
     1	   (1-AOMR(IW,IATP))   
         asmaer(IW)= asmaer(IW) + AGAR(IW,IATP)*NDS(IATP)
        END DO                       !aerosol type loop END

        albaer(IW)= albaer(IW)/NAER                             ! average back-scattering albedo
        asmaer(IW)= asmaer(IW)/Naer                             ! average asymmetry factor
	 
       END DO                        !wavelength loop END

      END IF                         !Naer ENDIF
       
      do IW=2,7
C---- Absorp & scattering at 550nm wavelength       
      if(wavei(iw-1).le.550.and.wavei(iw).ge.550) then
        tabsorp=absorp(iw-1)+(absorp(iw)-absorp(iw-1))/
     1     (wavei(iw)-wavei(iw-1))*(550-wavei(iw-1))
        tabsorp=tabsorp/1000.                            ! /km -> /m
        tscatter=extaer(iw-1)+(extaer(iw)-extaer(iw-1))/   ! scatter=extinction - absorption
     1    (wavei(iw)-wavei(iw-1))*(550-wavei(iw-1))-tabsorp
        tscatter=tscatter/1000.                          ! /km -> /m     
       return 
      endif
      enddo 
      END

      subroutine handle_err(iflag)
      print*,'error occurred in netcdf ',iflag
      stop
      end


       subroutine decomma(aline,pdata,ndata)  ! decode comma delimited data
       character*(*)aline
       real pdata(900)
       nlength=len_trim(aline)
       if(nlength.le.0) then
        ndata=0
	return
       endif 	

       ndata=1
       L2=index(aline,',')
       if(L2.eq.0) then   ! no comma
        read(aline,*)pdata(1)
        return
       else
        ndata=0	
       endif		
       L1=0
       do while(L2.lt.nlength)
	L2=index(aline(L1+1:nlength),',')+L1
	if(L2.eq.L1) L2=nlength+1
	ndata=ndata+1
	if(L2.eq.L1+1) then  ! empty one
	 pdata(ndata)=-9999.
	else
	 read(aline(L1+1:L2-1),*)pdata(ndata)
	endif 
	L1=L2
       enddo
       end

       subroutine decommac(aline,pcdata,ndata)  ! decode comma delimited data
       character*(*)aline
       character*80 pcdata(900)
       nlength=len_trim(aline)
       if(nlength.le.0) then
        ndata=0
	return
       endif 	

       ndata=1
       L2=index(aline,',')
       if(L2.eq.0) then   ! no comma
        pcdata(1)=aline
        return
       else
        ndata=0	
       endif		
       L1=0
       do while(L2.lt.nlength)
	L2=index(aline(L1+1:nlength),',')+L1
	if(L2.eq.L1) L2=nlength+1
	ndata=ndata+1
	if(L2.eq.L1+1) then  ! empty one
	 pcdata(ndata)=''
	else
	 pcdata(ndata)=aline(L1+1:L2-1)
	endif 
	L1=L2
       enddo
       end


!********************
! convert latitude longtitude to lambert conformal grid index in cross points
!

      subroutine mm5ll_xy(alat,alon,cenlat,cenlon,xp,yp,imm5,jmm5,dx)

	 real    alat,alon
         real    aa,an,dx,xoff,yoff,radeg,degn,degw,degp
         real    a1,r,xp,yp

         if(alon.lt.0) alon=alon+360
	 if(cenlon.lt.0) cenlon=cenlon+360


             aa    = 6370.997
             an    = 0.716

!            //
!            // xoff, and yoff are the grid distances from the center
!            // of the 12km SMM1 domain to the lower left corner of
!            // the domain for which emissions are being created.
!            //
             xoff  = (imm5 + 1) / 2.0 
             yoff  = (jmm5 + 1) / 2.0 

             radeg = 57.29577951

             c1    = abs(-cenlon - (90.0/an))
             c2    = aa/an*sin(30.0/radeg) *     
     1               (((tan((90.0-cenlat)/2.0/radeg)/tan(30.0/2.0/
     2                  radeg))**an))

             degn  = alat / radeg
             degw  = -alon / radeg
             degp  = -an * (degw + c1 / radeg)
             a1    = aa / dx / an / 2.0 / (tan(15.0 / radeg))**an
             r     = a1 * (tan((90.0 / radeg - degn) / 2.0))**an
             xp    = r * cos(degp) + xoff
             yp    = r * sin(degp) + c2 / dx + yoff

      return
      end

        subroutine ll2lc ( vals, grdlat, grdlon,  grdi, grdj) 

! Subroutine to convert from lat-lon to Lambert Conformal i,j.
! Provided by NRL Monterey; converted to C 6/15/94.
!                SUBROUTINE: ll2lc
!
!                PURPOSE: To compute i- and j-coordinates of a specified
!                         grid given the latitude and longitude points. 
!                         All latitudes in this routine start 
!                         with -90.0 at the south pole and increase 
!                         northward to +90.0 at the north pole.  The 
!			  longitudes start with 0.0 at the Greenwich 
!			  meridian and increase to the east, so that 
!			  90.0 refers to 90.0E, 180.0 is the inter- 
!			  national dateline and 270.0 is 90.0W. 
!
!		 INPUT VARIABLES: 
! 
!   vals+0	   reflat: latitude at reference point (iref,jref) 
!
!   vals+1	   reflon: longitude at reference point (iref,jref) 
! 
!   vals+2         iref:   i-coordinate value of reference point 
! 
!   vals+3	   jref:   j-coordinate value of reference point 
! 
!   vals+4	   stdlt1: standard latitude of grid 
! 
!   vals+5	   stdlt2: second standard latitude of grid (only required 
!		   if igrid = 2, lambert conformal) 
! 
!   vals+6	   stdlon: standard longitude of grid (longitude that 
!			   points to the north) 
! 
!   vals+7	   delx:   grid spacing of grid in x-direction 
!			   for igrid = 1,2,3 or 4, delx must be in meters 
!			   for igrid = 5, delx must be in degrees 
! 
!   vals+8	   dely:   grid spacing (in meters) of grid in y-direction 
!			   for igrid = 1,2,3 or 4, delx must be in meters 
!			   for igrid = 5, dely must be in degrees 
! 
!		   grdlat: latitude of point (grdi,grdj) 
! 
!		   grdlon: longitude of point (grdi,grdj)
! 
!		   grdi:   i-co ordinate(s) that this routine will generate 
!			   information for 
! 
!		   grdj:   j-coordinate(s) that this routine will generate 
!			   information for 

      real pi, pi2, pi4, d2r, r2d, radius, omega4
      real gcon,ogcon,ahem,deg,cn1,cn2,cn3,cn4,rih,xih,yih,rrih,check
      real alnfix,alon,x,y, vals(0:8),latref,lonref,iref,jref
   
      pi = 4.0*atan(1.0)
      pi2 = pi/2.0
      pi4 = pi/4.0
      d2r = pi/180.0
      r2d = 180.0/pi 
      radius = 6371229.0
      omega4 = 4.0*pi/86400.0

      latref = vals(0)
      lonref = vals(1)
      iref   = vals(2)    
      jref   = vals(3)    
      stdlt1 = vals(4)    
      stdlt2 = vals(5)    
      stdlon = vals(6)    
      delx   = vals(7)    
      dely   = vals(8)    
                    
! /* case where standard lats are the same */
! /* corrected by Dan Geiszler of NRL; fabs of the 
!   lats was required for shem cases */

      if(stdlt1.eq.stdlt2) then  
       gcon = sin(d2r*(abs(stdlt1)))
      else
       gcon = (log(sin((90.0-abs(stdlt1))*d2r)) 
     1   -log(sin((90.0-abs(stdlt2))*d2r)))
     2    /(log(tan((90.0-abs(stdlt1))*0.5*d2r))
     3     -log(tan((90.0-abs(stdlt2))*0.5*d2r)))
      endif
      ogcon = 1.0/gcon
      H = abs(stdlt1)/(stdlt1)              ! /* 1 for NHem, -1 for SHem */
      cn1 = sin((90.0-abs(stdlt1))*d2r)
      cn2 = radius*cn1*ogcon
      deg = (90.0-abs(stdlt1))*d2r*0.5
      cn3 = tan(deg)
      deg = (90.0-abs(latref))*d2r*0.5
      cn4 = tan(deg)
      rih = cn2*(cn4/cn3)**gcon

      xih =  rih*sin((lonref-stdlon)*d2r*gcon)
      yih = -rih*cos((lonref-stdlon)*d2r*gcon)*H
      deg = (90.0-grdlat*H)*0.5*d2r
      cn4 = tan(deg)
      rrih = cn2*(cn4/cn3)**gcon
      check  = 180.0-stdlon
      alnfix = stdlon+check
      alon   = grdlon+check

      if (alon<  0.0) alon = alon+360.0
      if (alon>360.0) alon = alon-360.0

      deg = (alon-alnfix)*gcon*d2r
      x =  rrih*sin(deg)
      y = -rrih*cos(deg)*H
      grdi = iref + (x-xih)/delx
      grdj = jref + (y-yih)/dely
      
      end

       subroutine decommach(aline,cdata,ndata)  ! decode space delimited character data
       character*(*)aline
       character*80 cdata(900)
             
       nlength=len_trim(aline)
       if(nlength.le.0) then
        ndata=0
        return
       endif    

       ndata=1
       L2=index(aline,',')
       if(L2.eq.0) then   ! no comma
        cdata(1)=trim(aline)
        return
       else
        ndata=0 
       endif            
       L1=0
       do while(L2.lt.nlength)
        L2=index(aline(L1+1:nlength),',')+L1
        if(L2.eq.L1) L2=nlength+1      ! near line end
        if(L2.gt.L1+1) then  ! non-empty one
         ndata=ndata+1
         cdata(ndata)=aline(L1+1:L2-1)
        endif 
 28     L1=L2
       enddo
       end
      subroutine lowcase(aline)
      character*(*) aline
      length=len_trim(aline)
      do n=1,length
       if(aline(n:n).le.'Z'.and.aline(n:n).ge.'A') then
        aline(n:n)=char(ichar(aline(n:n))-ichar('A')+ichar('a'))
c       else if(aline(n:n).eq.'-'.or.aline(n:n).eq.'['
c     1   .or.aline(n:n).eq.'+') then
c        aline(n:n)='_'
       endif
      enddo
      end
