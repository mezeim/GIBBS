c*****Cavity biased Gibbs ensemble Monte Carlo program
c     Written by Mihaly Mezei.
c     To dimension the program, execute the following changes:
c     change '499500' to '<(maxmol-1)*maxmol/2>' nol
c     change '1000' to '<maxmol>' nol
c     change '10000' to '<maxslv*maxmol>' nol
c     change '210' to '<maximumn number of rdf grids>' nol
c     change '10' to '<maxslv>' nol
c     change '1000000' to '<maximum number of cavities>' nol
c     change '100'  to '<maximum number of cavity grids in 1 dimension>' nol
      character*80 label
      common /id/ label(2),temp,fkt,cutof,cutof2,edgea,edge,hlfedg,
     -  prsag,prsar,pressg,pressr,bparg,bparr
      real*8 ug,ur,uug,uur,vg,vr,vvg,vvr,pg,pr,ppg,ppr,dva,dvag,dvar
      common /accu/ ug,ur,uug,uur,vg,vr,vvg,vvr,pg,pr,ppg,ppr,dva,
     -  dvag,dvar
      common /sreal/ eming,eminr,emaxg,emaxr,vol,volg,volr,vfacg,vfacr,
     -  cut2g,cut2r,vming,vmaxg,vminr,vmaxr
      common /spcpar/qqspc,qqspc2,qqspc4,alj,blj,alj6,alj66,blj12,blj122
      common /genslv/ c6(10),c12(10),cc6(10,10),cc12(10,10),qqg(10,10),
     -                siglj(10),epslj(10)
      common /svint/ nmc,nmaccg,nmaccr,neming,neminr,nemaxg,nemaxr,
     -  nmcmax,nmcrep,nrecd,nmcvch,nsync,niaccg,ndaccg,nmovg,nmovr,
     -  nvming,nvmaxg,npming,npmaxg,nvminr,nvmaxr,npminr,npmaxr,
     -  nvacc,nvaccg,nvaccr,nmccnt,ncnt
      common /info/ nslv,nmolg,nmolr,natg,natr,icopvg,icopvr,nmolmx
      common /grid/ ncover(100,100,100),igfree(1000000),ngfree,ngrid,
     -  ng1,ng2,ng3,fltngr,fltng3,grid,e2g,e2gm,rsppe2,rspme2
      common /coord/ relcrd(3,10),cslv(3,10),qslv(10),islv(10)
      common /dploy/ rsi(3,10),crmi(3),orj(3,3),cedgs,rotax,delvol
      common /ctotal/ cg(3,10000,2),cr(3,10000,2),org(3,3,1000),
     -  orr(3,3,1000),crmg(3,1000),crmr(3,1000)
      real*8 besg,besr,etotog,etotng,etotor,etotnr
      common /energy/ besg(1000,2),texg(1000,2),texog(1000,2),
     -                exg(499500,2),
     -                besr(1000,2),texr(1000,2),texor(1000,2),
     -                exr(499500,2),
     -                etotog,etotng,etotor,etotnr
      real*8 angcrg,angcrr
      common /distr/ angcrg(210),angcrr(210),nrg(210),nrr(210),
     -  rdbg(10,210),rdbr(10,210),angcbg(10,210),angcbr(10,210),
     -  ncnfbg(10),ncnfbr(10),rgrid,rgridi,rgrdig,rgrdir,rdmax2,rdmx2g,
     -  rdmx2r,rdm2ng,rdm2nr,rgring,rgrinr,nrgrid,icolcr,nconfg,nconfr,
     -  lenblk,nblk
      common /option/ iop(40)
      common /nmbrs/ autokc,autoa,rdtodg,boltz,pi,pi2
      character*4 iel
      common /iname/ iel(100)
      common /atwght/ avgd,aw(100)
      common /cavdat/ rsph,rsphs,corlng,edgrs,edgrs2,edgrq2
      common /cavprb/ ndistg(1000),ndistr(1000),pming,pmaxg,pminr,pmaxr
      common /cntrl/ ceg(100),cer(100),cvg(100),cng(100)
      common /seed/ ix
      common /texz/ tex0(1000)
      character*8 itg,itr
      common /label/ itg,itr
      dimension rno(4),iopn(40),islv1(10)
      real*8 ediff
      character*80 line
c     iop(1)  .eq. 1 : create inital configuration
c     iop(1)  .eq. 2 : input from unit 13, but update number of molecules
c     iop(2)  .eq. 0 : apply volume change to the r system
c     iop(2)  .eq. 1 : apply volume change to the g system
c     iop(3)  .eq. 0 : Gibbs ensemble simulation
c     iop(3)  .eq. 1 : two (T,P,N) ensemble simulations (iop(6)==>0)
c     iop(3)  .eq. 2 : two (T,V,mu) ensemble simulations
c     iop(4)  .eq. 0 : binary input on unit 13
c     iop(4)  .eq. 1 : alphanumeric input on unit 13
c     iop(6)  .eq. 0 : no particle exchange
c     iop(6)  .eq. 1 : use random insertion
c     iop(6)  .eq. 2 : use cavity-biased insertion
c     iop(8)  .eq. 0 : don't use the long-range correction
c     iop(8)  .eq. 1 : use the long-range correction
c     iop(9)  .eq. 0 : spherical cut-off for the electrostatic interactions
c     iop(9)  .eq. 1 : reaction field correction with eps(rf)=infinity   
c                      only if iop(11)=6
c     iop(10) .ne. 0 : restart, input from unit 10
c     iop(11) .eq. 0 : MCY water
c     iop(11) .eq. 1 : ST2 water
c     iop(11) .eq. 2 : SPC/TIPS water
c     iop(11) .eq. 3 : TIPS2/TIPS4P water
c     iop(11) .eq. 4 : LJ fluid
c     iop(11) .eq. 5 : general 12-6 (LJ) molecule
c     iop(11) .eq. 6 : general 12-6-1 molecule
c     iop(12) .eq. 0 : com synchronization at 10**5        -th MC step
c     iop(12) .gt. 0 : com synchronization at 10**(iop12-1)-th MC step
c     iop(13) .gt. 0 : echo inputted data
c     iop(15) .gt. 0 : save initial config on unit 13 and restart
c     iop(15) .eq. 2 : save initial config on unit 13 in alphanum and restart
c     iop(16) .gt. 0 : create two free-format coordinate file and stop
c     iop(17) .gt. 0 : calculate g(r) and angular correlations for the g sys
c     iop(17) .gt. 1 : calculate g(r) and angular correlations for the r sys
c     iop(17) .gt. 2 : calculate the absolute value of the correlations
c     The endings g and r in the variable names refer to the grid or
c     random insertion. grid insertion is done in the liquid and random
c     insertion in the vapor phase.
      character*18 dataname(9)
      data dataname /'options           ','                  ',
     -  'RDF info          ','system sizes      ',
     -  'step sizes        ','no of ats, run inf','atom data         ',
     -  'run length        ','restart run length'/
C@BM      call errset(208,256,-1,0,1,1)
      rdtodg=180.0/pi
      nmolmx=1000
      write (6,1041)
      call datprt
c     Read the options of the run
      irectyp=1
      read (5,1007,end=990) line
      read (line,1039,err=991) iop
      if(iop(9).gt.0.and.iop(11).ne.6) then
         write(6,*)' *** overwrite ***  no rf correction for uncharged',      
     *             'molecules!'
         write(6,*)   'iop(9) is changed to 0' 
         iop(9)=0
      endif
      if (iop(10). ne. 0) go to 8888
c     Read the verbal description of the run
      irectyp=2
      read (5,1000,end=990) label
      write (6,1040) label,iop
c     alj, blj: LJ parameters for SPC...
c     nrgrid: number of rdf girdpoints, rgrid: radial gridsize in ang
c     icolcr: the colum of the orientation matrix to use for correlation calc.
      irectyp=3
      read (5,1007,end=990) line
      read (line,2101,err=991) nrgrid,icolcr,rgrid,alj,blj
      if (iop(13) .gt. 0)
     -   write (6,1302) irectyp,nrgrid,icolcr,rgrid,alj,blj
      if (nrgrid .eq. 0) then
        nrgrid=210
        rgrid=0.05
      end if
      if (icolcr .eq. 0) icolcr=1
      rdmax2=(nrgrid*rgrid)**2
      rgridi=1.0/(rgrid)
c     Read the temperature in Kelvin and the cell size in a
c     coordinates are converted into a.u. for interval usage but volumes
c     are used in A**3 (as only their ratios are needed for updates)
      irectyp=4
      read (5,1007,end=990) line
      read (line,2103,err=991) nmmg,nmmr,temp,edgea,prsag,prsar
      if (iop(13) .gt. 0)
     -   write (6,1303) irectyp,nmmg,nmmr,temp,edgea,prsag,prsar
      pressg=prsag*avgd*1.e-27*(8.31441/0.0820575/4184.)
      pressr=prsar*avgd*1.e-27*(8.31441/0.0820575/4184.)
c     if nmmx<nmolx, the balance is deleted
c     if nmmx>nmolx, the balance of molecules is created randomly,
c     the new unit 13 is saved and the calculation stops.
      fkt=1.0/(boltz*temp)
      if (iop(11) .eq. 4) fkt=1.0/temp
c     Read the maximum translational and rotational displacements
      irectyp=5
      read (5,1007,end=990) line
      read (line,2101,err=991) i,i,cedgs,rotax,delvol
      if (iop(13) .gt. 0)
     -   write (6,1302) irectyp,i,i,cedgs,rotax,delvol
      rotax=rotax/rdtodg
c     Read the number of solvent atoms, number of cb grids, potential
c     cutoff (in a) and the hard core diameter
      irectyp=6
      read (5,1007,end=990) line
      read (line,2101,err=991) nslv,ngrid,cutof,rsph
      if (iop(13) .gt. 0)
     -   write (6,1302) irectyp,nslv,ngrid,cutof,rsph
      if (ngrid .eq. 0) ngrid=100
c     Atomic numbers and coordinates in a, local coord system
      qasum=0.0
      do i=1,nslv
        irectyp=7
        read (5,1007,end=990) line
        if (iop(11) .le. 4) then
          read (line,2102,err=991) islv(i),(cslv(j,i),j=1,3),qslv(i)
          if (iop(13) .gt. 0)
     -       write (6,1307) irectyp,islv(i),(cslv(j,i),j=1,3),qslv(i)
        else
          read (line,2102,err=991) 
     -      islv(i),(cslv(j,i),j=1,3),qslv(i),siglj(i),epslj(i)       
          if (iop(13) .gt. 0) write (6,1307)
     -       irectyp,islv(i),(cslv(j,i),j=1,3),qslv(i),siglj(i),epslj(i)
        end if
        qasum=qasum+abs(qslv(i))
        c6(i)=4.*epslj(i)*siglj(i)**6
        c12(i)=4.*epslj(i)*siglj(i)**12
      end do
      if (iop(11) .eq. 4 .and. nslv .gt. 1) then
        print *,'Lennard-Jones fluid can not have',nslv,' atoms'
        stop
      else if (iop(11) .eq. 6 .and. qasum .eq. 0.0) then
        print *,'Use type 5 for uncharged molecules' 
        if (iop(13) .eq. 0) then
          print *,'iop(11) is changed to 5'
          iop(11)=5
        end if
      end if
      if(iop(9).gt.0.and.iop(11).ne.6) then
         write(6,*)' *** Overwrite ***  no rf correction for uncharged',      
     *             'molecules!'
         write(6,*)   'iop(9) is changed to 0' 
         iop(9)=0
      endif
      call grelcd(cslv,relcrd,islv,nslv)
      call ewwint
      call init
      icopvg=1
      icopvr=1
      vfacg=1.0
      vfacr=1.0
      if (iop(1) .ne. 1) then
        if (iop(4) .eq. 0) read (13) nmolg,nmolr,ix,icopvv,volg,volr
        if (iop(4) .eq. 1) read (13,1002)nmolg,nmolr,ix,icopvv,volg,volr
        natg=nmolg*nslv
        natr=nmolr*nslv
        rewind 13
        if (iop(4) .eq. 0) then
          read (13) nmolg,nmolr,ix,icopvv,volg,volr,
     -      ((crmg(k,im),k=1,3),im=1,nmolg),
     -      ((crmr(k,im),k=1,3),im=1,nmolr)
          read (13,end=112)
     -      ((cg(k,ia,1),k=1,3),ia=1,natg),
     -      ((cr(k,ia,1),k=1,3),ia=1,natr)
        else
          read (13,1002) nmolg,nmolr,ix,icopvv,volg,volr
          read (13,1003)
     -      ((crmg(k,im),k=1,3),im=1,nmolg),
     -      ((crmr(k,im),k=1,3),im=1,nmolr)
          read (13,1003)
     -      ((cg(k,ia,1),k=1,3),ia=1,natg),
     -      ((cr(k,ia,1),k=1,3),ia=1,natr)
        end if
        go to 113
112     if (nslv .gt. 1) then
          write (6,8002) nmolg,nmolr,volg,volr
          stop
        end if
113     volo=(volg+volr)/2.0
        if (nslv .gt. 1) then
          vfacg=(volg/volo)**(1.0/3.0)
          vfacr=(volr/volo)**(1.0/3.0)
          e2=volo**(1.0/3.0)/2.0 + 0.0000001
          call findor(cg,crmg,org,nmolg,nslv,islv,relcrd,vfacg,e2,itg)
          call findor(cr,crmr,orr,nmolr,nslv,islv,relcrd,vfacr,e2,itr)
        end if
        if (abs(vol-volo)/(vol+volo) .gt. .0001 .and. iop(3).eq. 0) then
c         Transform to modified inital volumes
c         First scale the crm arrays
          vf=(vol/volo)**(1.0/3.0)
          do 110 im=1,nmolg
            do 110 k=1,3
110           crmg(k,im)=crmg(k,im)*vf
          do 111 im=1,nmolr
            do 111 k=1,3
111           crmr(k,im)=crmr(k,im)*vf
          if (iop(2) .eq. 1) then
c           Adjust the g system
            volgn=volg+2.0*(vol-volo)
            vfn=(volgn/vol)**(1.0/3.0)
            write (6,1031) vol,itg,volgn,volg,volr
            do 120 im=1,nmolg
              do 121 k=1,3
121             crmi(k)=crmg(k,im)*vfn
120           call putmol(im,cg,crmi,org)
            volg=volgn
            vfacg=(volg/vol)**(1.0/3.0)
            vfacr=(volr/vol)**(1.0/3.0)
          else
c           Adjust the r system
            volrn=volr+2.0*(vol-volo)
            vfn=(volrn/vol)**(1.0/3.0)
            write (6,1031) vol,itr,volrn,volg,volr
            do 122 im=1,nmolr
              do 123 k=1,3
123             crmi(k)=crmr(k,im)*vfn
122           call putmol(im,cr,crmi,orr)
            volr=volrn
            vfacg=(volg/vol)**(1.0/3.0)
            vfacr=(volr/vol)**(1.0/3.0)
          end if
        end if
        if (iop(1) .eq. 2) then
c         Add/subtract molecules based on nmmg and nmmr
          if (nmmg .gt. nmolg .and. nmmg .ne. 0)
     -      call genmol(nmolg,nmmg,nmolg,cg,crmg,org,edge,hlfedg,vfacg)
          if (nmmr .gt. nmolr .and. nmmr .ne. 0)
     -      call genmol(nmolr,nmmr,nmolr,cr,crmr,orr,edge,hlfedg,vfacr)
          if (nmmg .lt. nmolg .and. nmmg .ne. 0) then
            write (6,1005) nmolg,nmmg
            nmolg=nmmg
          end if
          if (nmmr .lt. nmolr .and. nmmr .ne. 0) then
            write (6,1005) nmolr,nmmr
            nmolr=nmmr
          end if
          natg=nmolg*nslv
          natr=nmolr*nslv
        end if
      else
c       Generate initial configuration
        volg=vol
        volr=vol
        call genmol(0,nmmg,nmolg,cg,crmg,org,edge,hlfedg,1.0)
        call genmol(0,nmmr,nmolr,cr,crmr,orr,edge,hlfedg,1.0)
        natg=nmolg*nslv
        natr=nmolr*nslv
        rewind 13
        if (iop(4) .eq. 0) then
          write (13) nmolg,nmolr,ix,icopvg,volg,volr,
     -  ((crmg(k,im),k=1,3),im=1,nmolg),((crmr(k,im),k=1,3),im=1,nmolr)
          write (13)
     -  ((cg(k,ia,1),k=1,3),ia=1,natg,1),((cr(k,ia,1),k=1,3),ia=1,natr)
        else
          write (13,1002) nmolg,nmolr,ix,icopvg,volg,volr
          write (13,1003) ((crmg(k,im),k=1,3),im=1,nmolg),
     -      ((crmr(k,im),k=1,3),im=1,nmolr)
          write (13,1003) ((cg(k,ia,1),k=1,3),ia=1,natg,1),
     -      ((cr(k,ia,1),k=1,3),ia=1,natr)
        end if
      end if
      irectyp=8
      read (5,1007,end=990) line
      read (line,1021,err=991) nmcmax,nmcrep,nrecd,nmcvch,nmccnt
      if (iop(13) .gt. 0)
     -   write (6,1308) irectyp,nmcmax,nmcrep,nrecd,nmcvch,nmccnt
      if (nmcvch .eq. 0) nmcvch=nmolg+nmolr
      lenblk=(nmcmax+1)/10
      if (lenblk .lt. nmcvch) lenblk=nmcvch
      if (iop(12) .eq. 0) then
        nsync=100000
      else
        nsync=10**(iop(12)-1)
      end if
      if (nmccnt .eq. 0) nmccnt=100000
100   call echopr
      call crorgn(edge)
      if (iop(6) .gt. 1)then 
        rsph=rsph/vfacg
        call initgr(crmg,nmolg,edge)
      end if
c     Initialize all accumulators
103   do 104 i=1,1000
        ndistg(i)=0
104     ndistr(i)=0
      if (iop(17) .gt. 0) then
        do 105 i=1,nrgrid
          nrg(i)=0
          nrr(i)=0
          angcrg(i)=0.d0
105       angcrr(i)=0.d0
        rdmx2g=rdmax2/vfacg**2
        rdmx2r=rdmax2/vfacr**2
        rgrdig=rgridi*vfacg
        rgrdir=rgridi*vfacr
      end if
      nmc=0
      cut2g=cutof2/vfacg**2
      cut2r=cutof2/vfacr**2
      nmaccg=0
      nmaccr=0
      nmovg=0
      nmovr=0
      niaccg=0
      ndaccg=0
      nvacc=0
      nvaccg=0
      nvaccr=0
      dva=0.d0
      dvag=0.d0
      dvar=0.d0
      ug=0.d0
      uug=0.d0
      vg=0.d0
      vvg=0.d0
      ng=0.d0
      nng=0.d0
      ur=0.d0
      uur=0.d0
      vr=0.d0
      vvr=0.d0
      nr=0.d0
      nnr=0.d0
      ncnt=0
      nblk=0
c     Calculation of  initial total energies
      call eism0(cg,crmg,exg,besg,etotog,nmolg,cut2g,vfacg,
     -  iop(17),nrg,angcrg,org,nconfg,icolcr,rgrdig,rdmx2g)
      call eism0(cr,crmr,exr,besr,etotor,nmolr,cut2r,vfacr,
     -  iop(17),nrr,angcrr,orr,nconfr,icolcr,rgrdir,rdmx2r)
      eming=etotog
      neming=nmc
      emaxg=etotog
      nemaxg=nmc
      eminr=etotor
      neminr=nmc
      emaxr=etotor
      nemaxr=nmc
      vming=volg
      vmaxg=volg
      nvming=nmc
      nvmaxg=nmc
      pming=nmolg
      pmaxg=nmolg
      npming=nmc
      npmaxg=nmc
      vminr=volr
      vmaxr=volr
      nvminr=nmc
      nvmaxr=nmc
      pminr=nmolr
      pmaxr=nmolr
      npminr=nmc
      npmaxr=nmc
      go to 8890
c     Restart
8888  call trnsfi(iopn,iop,40)
      call crstr(10)
      iop(13)=iopn(13)
      write (6,1040) label,iop
      if (iopn(16) .gt. 0) then
c       Create two free-format coordinate file
        do 8891 i=1,nslv
          islv1(i)=islv(i)
c         United atoms are printed as carbons
          if (islv1(i) .gt. 30) islv1(i)=6
8891    continue
        write (20,1001) ((iel(islv1(ia)),
     -    (cg(k,(im-1)*nslv+ia,icopvg),k=1,3),ia=1,nslv),im=1,nmolg)
        write (21,1001) ((iel(islv1(ia)),
     -    (cr(k,(im-1)*nslv+ia,icopvr),k=1,3),ia=1,nslv),im=1,nmolr)
        print *,'Configurations written on units 20 and 21 in free-form'
        stop
      end if
      if (iopn(15) .gt. 0) then
c       Save c on unit 13  and start from scratch
        natg=nmolg*nslv
        natr=nmolr*nslv
        rewind 13
        if (iopn(4) .eq. 0) then
          write (13) nmolg,nmolr,ix,icopvg,volg,volr,
     -      ((crmg(k,im),k=1,3),im=1,nmolg),
     -      ((crmr(k,im),k=1,3),im=1,nmolr)
          write (13)
     -      ((cg(k,ia,icopvg),k=1,3),ia=1,natg),
     -      ((cr(k,ia,icopvr),k=1,3),ia=1,natr)
        else
          write (13,1002) nmolg,nmolr,ix,icopvg,volg,volr
          write (13,1003)
     -      ((crmg(k,im),k=1,3),im=1,nmolg),
     -      ((crmr(k,im),k=1,3),im=1,nmolr)
          write (13,1003)
     -      ((cg(k,ia,icopvg),k=1,3),ia=1,natg),
     -      ((cr(k,ia,icopvr),k=1,3),ia=1,natr)
        end if
        write (6,8004) nmc,nmolg,nmolr,temp
        if (iopn(15) .eq. 2) write (6,8005)
        if (iop(3) .eq. 1) then
          edgeo=((volg+volr)/2.0)**(1.0/3.0)
          write (6,8003) edgeo
        end if
        stop
      end if
      irectyp=9
      read (5,1007,end=990) line
      read (line,1021,err=991) nmcmax,nmcrep
      if (iop(13) .gt. 0)
     -   write (6,1308) irectyp,nmcmax,nmcrep
      call ewwint
      call echopr
      write (6,8001) nmc
8890  nmolec=nmolg+nmolr
      volmin=(2.0*cutof)**3
      write (6,1029) etotog,etotor,nmolg,nmolr,nmolec,volg,volr,vol
      rnng=0.0
      rnnr=0.0
      if (nmolg .gt. 0) then
        npair=0
        nx=((nmolg-1)*nmolg)/2
        do 130 i=1,nx
          if (exg(i,icopvg) .ne. 0.0) npair=npair+1
130     continue
        rnng=float(2*npair)/float(nmolg)
      end if
      if (nmolr .gt. 0) then
        npair=0
        nx=((nmolr-1)*nmolr)/2
        do 131 i=1,nx
          if (exr(i,icopvr) .ne. 0.0) npair=npair+1
131     continue
        rnnr=float(2*npair)/float(nmolr)
      end if
      rnneg=nmolg*4.0*pi*cutof**3/3.0/volg
      rnner=nmolr*4.0*pi*cutof**3/3.0/volr
      write (6,1006) rnng,rnnr,rnneg,rnner
      if (nmcmax .gt. 0 .or. iop(13) .gt. 0) call stest(nmolg,nmolr,3)
      if (nmcmax .eq. 0 .and. iop(13) .gt. 0) call disfin
      if (nmcmax .eq. 0) stop
c     Begin Monte Carlo loop
      do 9999 iii=1,nmcmax
        nmc=nmc+1
        call random(2,rno)
        isg=nmolg*rno(1)+1
        isr=nmolr*rno(2)+1
        if (iop(13) .gt. 0 .and. mod(nmc,nmcrep) .eq. 0)
     -    call stest(isg,isr,iop(13))
        call move(isg,nmolg,cg(1,1,icopvg),crmg,org,etotog,etotng,
     -    iop(6),cut2g,ug,uug,eming,emaxg,neming,nemaxg,besg(1,icopvg),
     -    texg(1,icopvg),texog(1,icopvg),exg(1,icopvg),nmaccg,vfacg)
        call move(isr,nmolr,cr(1,1,icopvr),crmr,orr,etotor,etotnr,0,
     -    cut2r,ur,uur,eminr,emaxr,neminr,nemaxr,besr(1,icopvr),
     -    texr(1,icopvr),texor(1,icopvr),exr(1,icopvr),nmaccr,vfacr)
        if (mod(nmc,nmcvch) .eq. 0)  then
          if (iop(3) .eq. 0) then
c           Volume change attempt
            irej=-1
            call random(2,rno)
            dvol=delvol*(rno(1)-0.5)
            volng=volg+dvol
            volnr=2.0*vol-volng
            icpvng=3-icopvg
            icpvnr=3-icopvr
            call volcha(vol,volng,cg(1,1,icopvg),cg(1,1,icpvng),crmg,
     -        nmolg,nslv,vfacg,vfacng)
            call volcha(vol,volnr,cr(1,1,icopvr),cr(1,1,icpvnr),crmr,
     -        nmolr,nslv,vfacr,vfacnr)
            cut2ng=cutof2/vfacng**2
            cut2nr=cutof2/vfacnr**2
            if (iop(17) .gt. 0) then
              rdm2ng=rdmax2/vfacg**2
              rdm2nr=rdmax2/vfacr**2
              rgring=rgridi*vfacg
              rgrinr=rgridi*vfacr
            end if
            call eism0(cg(1,1,icpvng),crmg,exg(1,icpvng),besg(1,icpvng),
     -        etotng,nmolg,cut2ng,vfacng,
     -        iop(17),nrg,angcrg,org,nconfg,icolcr,rgring,rdm2ng)
            call eism0(cr(1,1,icpvnr),crmr,exr(1,icpvnr),besr(1,icpvnr),
     -        etotnr,nmolr,cut2nr,vfacnr,
     -        iop(17),nrr,angcrr,orr,nconfr,icolcr,rgrinr,rdm2nr)
            ediff=(etotog-etotng+etotor-etotnr)*fkt
            if (ediff .gt. 170.d0) then
              accpr=1.0
            else
              accpr=dexp(ediff+fkt*corlng*dvol*
     -          (nmolg**2/(volg*volng)-nmolr**2/(volr*volnr)))*
     -          (volng/volg)**nmolg*(volnr/volr)**nmolr
            end if
            if (accpr .gt. rno(2)) then
c             Volume change is accepted
              irej=+1
              icopvg=icpvng
              icopvr=icpvnr
              nvaccg=nvaccg+1
              nvaccr=nvaccr+1
              dva=dva+abs(dvol)
              etotog=etotng
              etotor=etotnr
              vfacg=vfacng
              vfacr=vfacnr
              volg=volng
              volr=volnr
              cut2g=cut2ng
              cut2r=cut2nr
              if (iop(17) .gt. 0) then
                rdmx2g=rdm2ng
                rdmx2r=rdm2nr
                rgrdig=rgring
                rgrdir=rgrinr
              end if
              call mnmxup(volg,nmc,vming,vmaxg,nvming,nvmaxg)
              call mnmxup(volr,nmc,vminr,vmaxr,nvminr,nvmaxr)
            end if
          else if (iop(3) .eq. 1) then
            irej=-1
            call random(4,rno)
            dvolg=delvol*(rno(1)-0.5)
            volng=volg+dvolg
            dvolr=delvol*(rno(2)-0.5)
            volnr=volr+dvolr
            icpvng=3-icopvg
            icpvnr=3-icopvr
            call volcha(vol,volng,cg(1,1,icopvg),cg(1,1,icpvng),crmg,
     -        nmolg,nslv,vfacg,vfacng)
            call volcha(vol,volnr,cr(1,1,icopvr),cr(1,1,icpvnr),crmr,
     -        nmolr,nslv,vfacr,vfacnr)
            cut2ng=cutof2/vfacng**2
            cut2nr=cutof2/vfacnr**2
            if (iop(17) .gt. 0) then
              rdm2ng=rdmax2/vfacg**2
              rdm2nr=rdmax2/vfacr**2
              rgring=rgridi*vfacg
              rgrinr=rgridi*vfacr
            end if
            call eism0(cg(1,1,icpvng),crmg,exg(1,icpvng),
     -        besg(1,icpvng),etotng,nmolg,cut2ng,vfacng,
     -        iop(17),nrg,angcrg,org,nconfg,icolcr,rgring,rdm2ng)
            call eism0(cr(1,1,icpvnr),crmr,exr(1,icpvnr),
     -        besr(1,icpvnr),etotnr,nmolr,cut2nr,vfacnr,
     -        iop(17),nrr,angcrr,orr,nconfr,icolcr,rgrinr,rdm2nr)
            ediff=(-pressg*dvolg+etotog-etotng)*fkt
            if (ediff .gt. 170.) then
              accpr=1.0
            else
              accpr=dexp(ediff)*(volng/volg)**nmolg
            end if
            if (accpr .gt. rno(3)) then
c             Accept for g
              irej=irej+2
              icopvg=icpvng
              nvaccg=nvaccg+1
              etotg=etotng
              vfacg=vfacng
              volg=volng
              cut2g=cut2ng
              if (iop(17) .gt. 0) then
                rdmx2g=rdm2ng
                rgrdig=rgring
              end if
              call mnmxup(volg,nmc,vming,vmaxg,nvming,nvmaxg)
              dvag=dvag+abs(dvolg)
            end if
            ediff=(-pressr*dvolr+etotor-etotnr)*fkt
            if (ediff .gt. 170.) then
              accpr=1.0
            else
              accpr=dexp(ediff)*(volnr/volr)**nmolr
            end if
            if (accpr .gt. rno(4)) then
c             Accept for r
              irej=irej+3
              icopvr=icpvnr
              nvaccr=nvaccr+1
              etotr=etotnr
              vfacr=vfacnr
              volr=volnr
              cut2r=cut2nr
              if (iop(17) .gt. 0) then
                rdmx2r=rdm2nr
                rgrdir=rgrinr
              end if
              call mnmxup(volr,nmc,vminr,vmaxr,nvminr,nvmaxr)
              dvar=dvar+abs(dvolr)
            end if
          end if
c         Accumulate
          vg=vg+volg
          vvg=vvg+volg*volg
          vr=vr+volr
          vvr=vvr+volr*volr
          ug=ug+etotog
          uug=uug+etotog**2
          ur=ur+etotor
          uur=uur+etotor**2
          if (mod(nmc,nmcrep) .eq. 0) then
            nmcv=nmc/nmcvch
            vgav=vg/nmcv
            vrav=vr/nmcv
            nvmn=nvming/1000
            nvmx=nvmaxg/1000
            nvacc=(nvaccg+nvaccr)/2
            pacc=float(nvacc)/float(nmcv)
            write (6,1023) nmc,etotng,etotnr,vgav,vrav,volng,volnr,
     -         vming,nvmn,vmaxg,nvmx,pacc,irej
          end if
          if (volng .lt. volmin .or. volnr .lt. volmin) then
            write (6,1030) volmin
            call csave(10)
            call disfin
            stop
          end if
        end if
        if (iop(6) .gt. 0) then
c         Particle exchange attempt
          indel=2*mod(nmc,2)-1
          if (indel .eq. 1 .and. nmolr .gt. 0 .or. 
     -        indel .eq. -1 .and. nmolg .gt. 0) then
            nmolng=nmolg
            call indlgr(cg(1,1,icopvg),crmg,org,texg(1,icopvg),
     -        texog(1,icopvg),besg(1,icopvg),etotog,etotng,nmolng,cut2g,
     -        indel,imolg,iop(6),pcavg,vfacg)
            nmolnr=nmolr
            call indlgr(cr(1,1,icopvr),crmr,orr,texr(1,icopvr),
     -        texor(1,icopvr),besr(1,icopvr),etotor,etotnr,nmolnr,cut2r,
     -        -indel,imolr,0,pcavr,vfacr)
            ediff=(etotog-etotng+etotor-etotnr)*fkt
            if (ediff .gt. 75.d0) then
              accpr=1.0
            else
              if (indel .eq. 1) then
                accpr=dexp(ediff+((2*nmolr-1)*corlng/volr-(2*nmolg-1)*
     -            corlng/volg)*fkt)*(volg*nmolr)/(volr*nmolng)*pcavg
              else
                accpr=dexp(ediff+((2*nmolg-1)*corlng/volg-(2*nmolr-1)*
     -            corlng/volr)*fkt)*(volr*nmolg)/(volg*nmolnr)/pcavg
              end if
            end if
          else
            accpr=0.0
          end if
          call random(1,rno)
          if (accpr .gt. rno(1)) then
c           Particle exchange is accepted
            etotog=etotng
            etotor=etotnr
            nmolg=nmolng
            nmolr=nmolnr
            p=nmolg
            call mnmxup(p,nmc,pming,pmaxg,npming,npmaxg)
            p=nmolr
            call mnmxup(p,nmc,pminr,pmaxr,npminr,npmaxr)
            if (indel .eq. 1) then
c             Insert into g, delete from r
              niaccg=niaccg+1
              if (iop(6) .gt. 1) call addcov(crmg(1,imolg))
              call savegy(besg(1,icopvg),texg(1,icopvg),tex0,
     -          exg(1,icopvg),nmolg,nmolg)
              nmolr1=nmolr+1
c             Deduct deleted molec-s contributions from bes
              call getegy(texor(1,icopvr),exr(1,icopvr),imolr,nmolr1)
              do 500 i=1,nmolr1
500             besr(i,icopvr)=besr(i,icopvr)-texor(i,icopvr)
              if (imolr .le. nmolr) then
                call movex(texor(1,icopvr),exr(1,icopvr),imolr,nmolr1)
                besr(imolr,icopvr)=besr(nmolr1,icopvr)
                call trnsfr(crmr(1,imolr),crmr(1,nmolr1),3)
                if (nslv .gt. 1) then
                  call trnsfr(cr(1,(imolr-1)*nslv+1,icopvr),
     -              cr(1,nmolr*nslv+1,icopvr),3*nslv)
                  call trnsfr(orr(1,1,imolr),orr(1,1,nmolr1),9)
                end if
              end if
            else
c             Insert into r, remove for g
              ndaccg=ndaccg+1
              call savegy(besr(1,icopvr),texr(1,icopvr),tex0,
     -          exr(1,icopvr),nmolr,nmolr)
              nmolg1=nmolg+1
c             Deduct deleted molec-s contributions from bes
              call getegy(texog(1,icopvg),exg(1,icopvg),imolg,nmolg1)
              do 510 i=1,nmolg1
510             besg(i,icopvg)=besg(i,icopvg)-texog(i,icopvg)
              if (imolg .le. nmolg) then
                call movex(texog(1,icopvg),exg(1,icopvg),imolg,nmolg1)
                besg(imolg,icopvg)=besg(nmolg1,icopvg)
                call trnsfr(crmg(1,imolg),crmg(1,nmolg1),3)
                if (nslv .gt. 1) then
                  call trnsfr(cg(1,(imolg-1)*nslv+1,icopvg),
     -              cg(1,nmolg*nslv+1,icopvg),3*nslv)
                  call trnsfr(org(1,1,imolg),org(1,1,nmolg1),9)
                end if
              end if
            end if
          else
c           Exchange rejected
            if (indel .eq. -1 .and. iop(6) .gt. 1)
     -        call addcov(crmg(1,imolg))
            imolg=-imolg
            imolr=-imolr
          end if
c         Accumulate
          pg=pg+nmolg
          ppg=ppg+nmolg*nmolg
          pr=pr+nmolr
          ppr=ppr+nmolr*nmolr
          ug=ug+etotog
          uug=uug+etotog**2
          ur=ur+etotor
          uur=uur+etotor**2
          if (nmolg .gt. 0) ndistg(nmolg)=ndistg(nmolg)+1
          if (nmolr .gt. 0) ndistr(nmolr)=ndistr(nmolr)+1
          if (mod(nmc,nmcrep) .eq. 0) then
            pgav=pg/nmc
            prav=pr/nmc
            nmn=pming
            nmx=pmaxg
            nnmn=npming/1000
            nnmx=npmaxg/1000
            pacc=float(niaccg+ndaccg)/float(nmc)
            write (6,1022) nmc,etotng,etotnr,pgav,prav,nmolng,nmolnr,
     -         nmn,nnmn,nmx,nnmx,pacc,imolg,imolr
          end if
        end if
        if (mod(nmc,nmccnt) .eq. 0 .and. ncnt .lt. 100) then
c         Save control function info on energy, volume, nmol
          ncnt=ncnt+1
          ceg(ncnt)=ug
          cvg(ncnt)=vg
          cng(ncnt)=pg
          cer(ncnt)=ur
        end if
        if (iop(17) .gt. 0) then
          if (mod(nmc,lenblk) .eq. 0) then
c           Save control function info on rdf, corr 
            if (nblk .eq. 10) then
c             Condense information
              lenblk=lenblk*2
              print *,'RDF control function data was condensed. new',
     -          ' block length=',lenblk,' MC steps'
              nblk=nblk/2
              do 550 ib=1,5
                ncnfbg(ib)=ncnfbg(2*ib)
                ncnfbr(ib)=ncnfbr(2*ib)
                do 550 ig=1,210
                  rdbg(ib,ig)=rdbg(2*ib,ig)
                  rdbr(ib,ig)=rdbr(2*ib,ig)
                  angcbg(ib,ig)=angcbg(2*ib,ig)
550               angcbr(ib,ig)=angcbr(2*ib,ig)
            else
c             Save information
              nblk=nblk+1
              ncnfbg(nblk)=nconfg
              ncnfbr(nblk)=nconfr
              do 551 ig=1,210
                rdbg(nblk,ig)=nrg(ig)
                rdbr(nblk,ig)=nrr(ig)
                angcbg(nblk,ig)=angcrg(ig)
551             angcbr(nblk,ig)=angcrr(ig)
            end if
          end if
        end if
        if (mod(nmc,nrecd) .eq. 0 .or. iii .eq. nmcmax) call csave(10)
        if (mod(nmc,nsync) .eq. 0) then
c         Synchronize com with atomic coordinates
          devmxg=0.0
          if (nmolg .gt. 0) then
            do 600 im=1,nmolg
              call cofms(cg(1,(im-1)*nslv+1,icopvg),crmi,islv,1,nslv)
              do 604 k=1,3
604             crmi(k)=crmi(k)/vfacg
              dev2=0.0
              do 602 k=1,3
602             dev2=dev2+(crmi(k)-crmg(k,im))**2
              if (dev2 .gt. devmxg) devmxg=dev2
              do 601 k=1,3
601             crmi(k)=crmg(k,im)*vfacg
600           call putmol(im,cg,crmi,org(1,1,icopvg))
            call eism0(cg(1,1,icopvg),crmg,exg(1,icopvg),
     -        besg(1,icopvg),etotog,nmolg,cut2g,vfacg,
     -        0,nrg,angcrg,org,nconfg,icolcr,rgrdig,rdmx2g)
          end if
          devmxr=0.0
          if (nmolr .gt. 0) then
            do 610 im=1,nmolr
              call cofms(cr(1,(im-1)*nslv+1,icopvr),crmi,islv,1,nslv)
              do 614 k=1,3
614             crmi(k)=crmi(k)/vfacr
              dev2=0.0
              do 603 k=1,3
603             dev2=dev2+(crmi(k)-crmr(k,im))**2
              if (dev2 .gt. devmxr) devmxr=dev2
              do 611 k=1,3
611             crmi(k)=crmr(k,im)*vfacr
610           call putmol(im,cr,crmi,orr(1,1,icopvr))
            call eism0(cr(1,1,icopvr),crmr,exr(1,icopvr),
     -        besr(1,icopvr),etotor,nmolr,cut2r,vfacr,
     -        0,nrg,angcrr,orr,nocnfr,icolcr,rgrdir,rdmx2r)
          end if
          devmxg = sqrt(devmxg)
          devmxr = sqrt(devmxr)
          write (6,1009) devmxg,devmxr
        end if
        if (mod(nmc,1000000) .eq. 0 .and. iii .lt. nmcmax) call disfin
9999  continue
      if (iop(13) .eq. 0 .and. nmc .gt. 1000) call stest(isg,isr,3)
      call disfin
      stop
990   write (6,1011) dataname(irectyp)
      stop
991   write (6,1012) dataname(irectyp),line
      stop
607   format(/,' First starting point = ',i7,/)
1000  format(a80)
1001  format(1x,a4,2x,3f15.5) 
1002  format(2i6,i12,i2,2e15.7)
1003  format(3f20.10)
1005  format(' Number of molecules=',i4,' reduced to ',i4,' as asked')
1006  format(' Average number for neighbours within cutoff distance',
     -  ' for the g and r systems=',2f8.2,' estimated:',2f8.2)
1009  format(' C.O.M. synchronized with atoms - max. dev=',2f12.6)
1011  format(' ***** ERROR: run out of data while trying to read ',a)
1012  format(' ***** ERROR: invalid input for ',a,' the line read:',/,a)
1007  format(a)
1021  format(7i10)
1022  format(' N=',i7,' e=',e12.6,1x,e12.6,' <n>=',f6.2,1x,f6.2,
     -  ' ng,r=',i4,1x,i4,' nmn=',i4,
     -  '(',i4,'k) nmx=',i4,'(',i4,'k) a=',f7.5,' is=',2i4)
1023  format(' N=',i7,' e=',e12.6,1x,e12.6,' <v>=',f8.2,1x,f8.1,
     -  ' vg,r=',f8.2,1x,f8.1,' vmn=',f6.0,
     -  '(',i4,'k) vmx=',f6.0,'(',i4,'k) a=',f5.3,
     -  ' r=',i2)
1029  format(/,' Total energies at the starting point =',2e14.6,/,
     -  ' number of molecules=',3i5,/,' current volumes=',2f15.5,
     -  ' initial common volume=',f12.5,' A**3')
1030  format(' Volume dropped below the minimum volume as defined by',
     -   ' the cutoff sphere (',f10.3,' A**3)')
1031  format(' Initial common volume changed to ',f10.3,' and the ',
     -  a8,' volume changed to ',f10.3,' A**3',/,
     -  8x,'Previous volumes=',2f12.3)
1039  format(40i2)
1041  format(1h1,//' Gibbs-ensemble Monte Carlo program',
     -  ' using the cavity-biased algorithm',/,
     -  ' written by Mihaly Mezei',//,
     -  ' Maximum number of molecules=1000',/,
     -  ' Maximum number of atoms/molecule=10',/,
     -  ' Maximum number of cavity grids=100**3',/,
     -  ' Maximum number of cavities=1000000',/)
1040  format(/,2(1x,a80,/),/,' Options:',40i3,/)
2101  format(2i5,3f20.0)
2102  format(i5,3f15.10,3f10.0)
2103  format(2i5,4f15.5)
8001  format(/,' Restart run from MC step no ',i7,/)
8002  format(1x,20('*'),' ERROR: unit 13 lacks coordinate array record',
     -  '  nmolg,nmolr=',2i6,' volg,volr=',2f12.5)
8003  format(' use ',f12.6,' A edge for Gibbs run')
8004  format(/,' starting coordinates are saved on unit 13 ',
     -  'and accumulators are reset at nMC=',i7,' nmol=',2i5,
     -  ' temp=',f6.1,' Kelvin')
8005  format(' data saved in alphanumeric format')
1302  format(' RECTYP ',i2,':',2i5,3f20.10)
1303  format(' RECTYP ',i2,':',2i5,4f15.6)
1307  format(' RECTYP ',i2,':',i5,3f15.10,3f10.4)
1308  format(' RECTYP ',i2,':',7i10)
      end
      block data
      common /nmbrs/ autokc,autoa,rdtodg,boltz,pi,pi2
      character*4 iel
      common /iname/ iel(100)
      common /atwght/ avgd,aw(100)
      common /seed/ ix
      common /texz/ tex0(1000)
      character*8 itg,itr
      common /label/ itg,itr
      data pi/3.141593/,pi2/6.2831853/
      data autokc/627.55/,boltz/.001987/,autoa/0.529167/,avgd/6.023e+23/
c     Parameters for the water-water potential
      data iel/'   H','  He','  Li','  Be','   B','   C','   N',
     -  '   O','   F','  Ne','  Na','  Mg','  Al','  Si','   P',
     -  '   S','  Cl','  Ar','   K',11*'    ','CH1 ','CH2 ','CH3 ',
     -  65*'    ','  LJ','  LP'/
      data aw/1.00797,4.0026,6.939,9.0122,10.811,12.01115,
     -  14.0067,15.9994,18.9984,20.17,22.98977,24.305,
     -  26.98154,28.08,30.97376,32.06,35.453,39.94,39.09,11*0.0,
     -  13.01912,14.02709,15.03506,65*0.0,1.0,0.0/
      common /cnew/ qq,dqq,a1,a2,a3,a4,b1,b2,b3,b4,ab1,ab2,ab3,ab4
      data qq/0.514783/,dqq/1.029566/,a1/1734.196/,a2/1.061887/,
     -     a3/2.319395/,a4/0.436006/,b1/2.726696/,b2/1.460975/,
     -     b3/1.567367/,b4/1.181792/
c*****Random number seed initialization
      data ix/1357/
      data tex0/1000*0.0/
      data itg/'G-system'/,itr/'R-system'/
      end
      subroutine stest(isg,isr,iop13)
c*****Self test on stored energies
      character*80 label
      common /id/ label(2),temp,fkt,cutof,cutof2,edgea,edge,hlfedg,
     -  prsag,prsar,pressg,pressr,bparg,bparr
      common /sreal/ eming,eminr,emaxg,emaxr,vol,volg,volr,vfacg,vfacr,
     -  cut2g,cut2r,vming,vmaxg,vminr,vmaxr
      common /svint/ nmc,nmaccg,nmaccr,neming,neminr,nemaxg,nemaxr,
     -  nmcmax,nmcrep,nrecd,nmcvch,nsync,niaccg,ndaccg,nmovg,nmovr,
     -  nvming,nvmaxg,npming,npmaxg,nvminr,nvmaxr,npminr,npmaxr,
     -  nvacc,nvaccg,nvaccr,nmccnt,ncnt
      common /coord/ relcrd(3,10),cslv(3,10),qslv(10),islv(10)
      common /ctotal/ cg(3,10000,2),cr(3,10000,2),org(3,3,1000),
     -  orr(3,3,1000),crmg(3,1000),crmr(3,1000)
      real*8 besg,besr,etotog,etotng,etotor,etotnr
      common /energy/ besg(1000,2),texg(1000,2),texog(1000,2),
     -                exg(499500,2),
     -                besr(1000,2),texr(1000,2),texor(1000,2),
     -                exr(499500,2),
     -                etotog,etotng,etotor,etotnr
      real*8 angcrg,angcrr
      common /distr/ angcrg(210),angcrr(210),nrg(210),nrr(210),
     -  rdbg(10,210),rdbr(10,210),angcbg(10,210),angcbr(10,210),
     -  ncnfbg(10),ncnfbr(10),rgrid,rgridi,rgrdig,rgrdir,rdmax2,rdmx2g,
     -  rdmx2r,rdm2ng,rdm2nr,rgring,rgrinr,nrgrid,icolcr,nconfg,nconfr,
     -  lenblk,nblk
      common /info/ nslv,nmolg,nmolr,natg,natr,icopvg,icopvr,nmolmx
      common /grid/ ncover(100,100,100),igfree(1000000),ngfree,ngrid,
     -  ng1,ng2,ng3,fltngr,fltng3,grid,e2g,e2gm,rsppe2,rspme2
      common /nmbrs/ autokc,autoa,rdtodg,boltz,pi,pi2
      common /option/ iop(40)
      character*8 itg,itr
      common /label/ itg,itr
      real*8 esmg,esmr,esmtg,esmtr,esmxg,esmxr
      dimension cm(3),z(3),cc(3),ortt(3,3)
      equivalence (iop11,iop(11))
      call csave(10)
      write (6,1111) nmc,isg,isr
c     Collect eij sum, find nearest distance
      esmg=0.d0
      nij=nmolg*(nmolg-1)/2
      do 101 i=1,nij
101     esmg=esmg+exg(i,icopvg)
      esmr=0.d0
      nij=nmolr*(nmolr-1)/2
      do 111 i=1,nij
111     esmr=esmr+exr(i,icopvr)
c     Collect tesi sum.
      esmtg=0.d0
      do 102 i=1,nmolg
102     esmtg=esmtg+besg(i,icopvg)
      esmtr=0.d0
      do 112 i=1,nmolr
112     esmtr=esmtr+besr(i,icopvr)
      esmtg=esmtg/2.0
      esmtr=esmtr/2.0
      write (6,1000) etotog,etotor,esmg,esmr,esmtg,esmtr
      call eism1(cg(1,1,icopvg),crmg,texg(1,icopvg),besg(1,icopvg),
     -  etotng,etotog,isg,nmolg,cut2g,vfacg)
      call eism1(cr(1,1,icopvr),crmr,texr(1,icopvr),besr(1,icopvr),
     -  etotnr,etotor,isr,nmolr,cut2r,vfacr)
      call getegy(texog,exg(1,icopvg),isg,nmolg)
      call getegy(texor,exr(1,icopvr),isr,nmolr)
      esmg=0.d0
      esmxg=0.d0
      do 103 i=1,nmolg
        if (i .ne. isg) esmg=esmg+texg(i,icopvg)
        if (i .ne. isg) esmxg=esmxg+texog(i,1)
103   continue
      esmr=0.d0
      esmxr=0.d0
      do 113 i=1,nmolr
        if (i .ne. isr) esmr=esmr+texr(i,icopvr)
        if (i .ne. isr) esmxr=esmxr+texor(i,1)
113   continue
      write (6,1004) isg,besg(isg,icopvg),isr,besr(isr,icopvr),
     -  esmxg,esmxr,esmg,esmr
      if (iop13 .gt. 2) then
        call eism0(cg(1,1,icopvg),crmg,exg(1,icopvg),besg(1,icopvg),
     -    etotng,nmolg,cut2g,vfacg,
     -    0,nrg,angcrg,org,nconfg,icolcr,rgrdig,rdmx2g)
        call eism0(cr(1,1,icopvr),crmr,exr(1,icopvr),besr(1,icopvr),
     -    etotnr,nmolr,cut2r,vfacr,
     -    0,nrr,angcrr,orr,nconfr,icolcr,rgrdir,rdmx2r)
        rming=0.0
        iming=0
        jming=0
        if (nmolg .gt. 1) then
          rming=1.e+20
          do 120 i=2,nmolg
            j1=i-1
            do 120 j=1,j1
              do 121 k=1,3
121             z(k)=crmg(k,i)-crmg(k,j)
              call cellfn(cc(1),cc(2),cc(3),z(1),z(2),z(3),hlfedg,ipbc)
              if (ipbc .gt. 1) then
                do 122 k=1,3
122               z(k)=z(k)+cc(k)
              end if
              r2=0.0
              do 123 k=1,3
123             r2=r2+z(k)**2
              if (r2 .lt. rming) then
                rming=r2
                iming=i
                jming=j
              end if
120       continue
          rming=sqrt(rming)*vfacg
        end if
        rminr=0.0
        iminr=0
        jminr=0
        if (nmolr .gt. 1) then
          rminr=1.e+20
          do 130 i=2,nmolr
            j1=i-1
            do 130 j=1,j1
              do 131 k=1,3
131             z(k)=crmr(k,i)-crmr(k,j)
              call cellfn(cc(1),cc(2),cc(3),z(1),z(2),z(3),hlfedg,ipbc)
              if (ipbc .gt. 1) then
                do 132 k=1,3
132               z(k)=z(k)+cc(k)
              end if
              r2=0.0
              do 133 k=1,3
133             r2=r2+z(k)**2
              if (r2 .lt. rminr) then
                rminr=r2
                iminr=i
                jminr=j
              end if
130       continue
          rminr=sqrt(rminr)*vfacr
        end if
        write (6,1002) etotng,etotnr,rming,iming,jming,rminr,iminr,jminr
        call crstr(10)
      end if
      if (iop(6) .lt. 2 .or. iop13 .lt. 2) go to 300
c     Check cover list for consistency
      ncerr=0
      do 61 ind=1,ngfree
        indx=igfree(ind)
        ix1=indx-1
        indx1=ix1/ng1
        iz=indx-indx1*ng1
        ix=indx1/ng1
        iy=indx1-ix*ng1+1
        ix=ix+1
        if (ncover(ix,iy,iz) .eq. -ind) go to 61
        ncerr=ncerr+1
        write (6,1009) ind,igfree(ind),ix,iy,iz,ncover(ix,iy,iz)
61    continue
      write (6,1010) ncerr
      ncerr=0
      do 62 ic=1,ngrid
        do 62 jc=1,ngrid
          do 62 kc=1,ngrid
            if (ncover(ic,jc,kc) .gt. 0) go to 62
            indx=igfree(-ncover(ic,jc,kc))
            ix1=indx-1
            indx1=ix1/ng1
            iz=indx-indx1*ng1
            ix=indx1/ng1
            iy=indx1-ix*ng1+1
            ix=ix+1
            if (ic .eq. ix .and. jc .eq. iy .and. kc .eq. iz) go to 62
            ncerr=ncerr+1
            write (6,1011) ic,jc,kc,ncover(ic,jc,kc),ix,iy,iz
62    continue
      write (6,1012) ncerr
      if (ng3 .le. 1000000) then
        ngf0=ngfree
        do 65 i=1,nmolg
65        call remcov(crmg(1,i))
        write (6,1013) ngf0,ngfree
      end if
      call crstr(10)
c     Check for pbc errors
300   npbcer=0
      do 301 i=1,nmolg
        call cellfn(cc1,cc2,cc3,crmg(1,i),crmg(2,i),crmg(3,i),hlfedg,
     -    ipbc)
        if (ipbc .eq. 1) go to 301
        write (6,1005) i,ipbc,(crmg(k,i),k=1,3)
        npbcer=npbcer+1
301   continue
      ncmerr=0
c     Check c.o.m. vs the coordinates of atoms
      if (nmolg .gt. 0 .and. iop(11) .ne. 4) then
        do 310 im=1,nmolg
          call cofms(cg(1,(im-1)*nslv+1,icopvg),cm,islv,1,nslv)
          d=0.0
          do 311 k=1,3
311         d=d+abs(crmg(k,im)*vfacg-cm(k))
          if (d .gt. 0.001) then
            ncmerr=ncmerr+1
            write (6,1006) im,d,itg
          end if
310     continue
      end if
      if (nmolr .gt. 0 .and. iop(11) .ne. 4) then
        do 320 im=1,nmolr
          call cofms(cr(1,(im-1)*nslv+1,icopvr),cm,islv,1,nslv)
          d=0.0
          do 321 k=1,3
321         d=d+abs(crmr(k,im)*vfacr-cm(k))
          if (d .gt. 0.001) then
            ncmerr=ncmerr+1
            write (6,1006) im,d,itr
          end if
320     continue
      end if
      noerr=0
      if (nmolg. gt. 0 .and. iop(11) .ne. 4) then
        do 850 im=1,nmolg
            call getort(cg(1,(im-1)*nslv+1,icopvg),crmg(1,im),islv,
     -        relcrd,nslv,ortt)
          devs=0.0
          do 851 k=1,3
            do 851 l=1,3
851           devs=devs+abs(ortt(k,l)-org(k,l,im))
          if (devs .gt. 0.001) then
            noerr=noerr+1
            write (6,1007) im,devs,itg
          end if
850     continue
      end if
      if (nmolr. gt. 0 .and. iop(11) .ne. 4) then
        do 860 im=1,nmolr
            call getort(cr(1,(im-1)*nslv+1,icopvr),crmr(1,im),islv,
     -        relcrd,nslv,ortt)
          devs=0.0
          do 861 k=1,3
            do 861 l=1,3
861           devs=devs+abs(ortt(k,l)-orr(k,l,im))
          if (devs .gt. 0.001) then
            noerr=noerr+1
            write (6,1007) im,devs,itr
          end if
860     continue
      end if
      write (6,1003) npbcer,ncmerr,noerr
c     Random number generator test
      esmo=0.d0
      nless=0
      nerr=0
      call random(1000,texg)
      do 104 i=1,1000
        esmo=esmo+texg(i,1)
        if (texg(i,1) .lt. 0.5) nless=nless+1
        if (texg(i,1) .lt. 0.0 .or. texg(i,1) .gt. 1.0) nerr=nerr+1
104   continue
      nmore=1000-nless
      ravg=esmo/1000.d0
      write (6,1014) nless,nmore,ravg,nerr
900   call crstr(10)
      return
1000  format(' etoto(gr)=',2e18.9,/' eij sums =',2e18.9,/,
     -       ' tesi sums=',2e18.9)
1001  format(' tesi(',i4,')=',2e18.9,/,' recomputed=',2e18.9)
1002  format(' Completely recomputed etoto(g,r)=',2e18.9,/,
     -  ' closest approach (g)=',f9.3,' (',2i4,')  (r)=',f9.3,' (',2i4,
     -  ')')
1003  format(' The number of PBC, c.o.m. and ormat errors=',3i4)
1004  format(' besg(',i3,')=',e18.9,' besr(',i3,')=',e18.9,/,
     -       '  exg sum =',e18.9,'  exr sum =',e18.9,/,
     -       ' texg sum =',e18.9,' texr sum =',e18.9)
1006  format(' C.O.M. for molecule',i4,' is off by',f12.6,' (',a8,')')
1007  format(' Orient matrix (mol',i4,') is off by',f9.6,' (',a8,')')
1005  format(' PBC error: molecule',i4,' is in cell',i3,3f10.4)
1009  format(' ERROR: igfree(',i5,')=',i8,' points to ncover(',
     -  i3,',',i3,',',i3,')=',i8)
1010  format(' The number of igfree errors=',i6)
1011  format(' ERROR: ncover(',i3,',',i3,',',i3,')=',i8,' points back',
     -  ' to ',3i4)
1012  format(' The number of ncover errors=',i6)
1013  format(' The number of free grids=',i8,
     -  '   the number of free grids after removing all atoms=',i8)
1014  format(' Single random test, below and above 0.5:',2i5,' average='
     -  ,f9.6,' number of errors=',i5)
1015  format(' Double random test, below and above 0.5:',2i5,
     -  ' average=',f9.6,' number of errors=',i5)
1111  format(' Self test at nmc=',i8,' molecules ',2i4)
      end
      subroutine init
c*****Initialize constants of the run
      common /cavdat/ rsph,rsphs,corlng,edgrs,edgrs2,edgrq2
      common /sreal/ eming,eminr,emaxg,emaxr,vol,volg,volr,vfacg,vfacr,
     -  cut2g,cut2r,vming,vmaxg,vminr,vmaxr
      character*80 label
      common /id/ label(2),temp,fkt,cutof,cutof2,edgea,edge,hlfedg,
     -  prsag,prsar,pressg,pressr,bparg,bparr
      common /coord/ relcrd(3,10),cslv(3,10),qslv(10),islv(10)
      common /genslv/ c6(10),c12(10),cc6(10,10),cc12(10,10),qqg(10,10),
     -                siglj(10),epslj(10)
      common /option/ iop(40)
      common /info/ nslv,nmolg,nmolr,natg,natr,icopvg,icopvr,nmolmx
      common /nmbrs/ autokc,autoa,rdtodg,boltz,pi,pi2
      edge=edgea
c     Simple cubic
      vol=edgea**3
      hlfedg=edge/2.0
      if (cutof .eq. 0.0 .or. cutof .gt. edgea/2.0) cutof=edgea/2.0
c     Constants to generate random pts away from the wall
      edgrs2=hlfedg-rsph
      edgrs=edgrs2*2.0
200   ct=cutof
      cutof2=ct**2
      corlng=0.0
      if (iop(8) .gt. 0 .and. iop(11) .eq. 4) then
        corlng=8.0*pi*(ct**(-9)/9.0-ct**(-3)/3.0)
c       For 2-dimensional LJ:
c       corlng=4.0*pi*(ct**(-10)/10.0-ct**(-5)/5.0)
      else if (iop(8) .gt. 0 .and. iop(11) .gt. 4) then
c       Molecular LJ correction
        corlng=-elr(6,ct,1.0,nslv,nslv,relcrd,relcrd,cc6,10) +
     -         elr(12,ct,1.0,nslv,nslv,relcrd,relcrd,cc12,10) 
        qsum=0.0
        do 10 i=1,nslv
10        qsum=qsum+abs(qslv(i))
        if (qsum .gt. 0.001) print *,'----- WARNING: ',
     -     'the long-range correction will not include the charges'
      else if (iop(8) .gt. 0) then
        print *,
     -    'No long-range correction is available for this potential'
      end if
      return
      end
      function elr(k2,rc,den,n1,n2,d1,d2,cij,ndim)
c     Calculate the long-range contribution from 1/r**k2 site-site interactions
c     between molecules d1 and d2, d1 being the central. den is the density of
c     d2. 
      dimension d1(3,n1),d2(3,n2),cij(ndim,ndim)
      dimension dd1(10),dd2(10)
      common /nmbrs/ autokc,autoa,rdtodg,boltz,pi,pi2
      elr=0.0
      k=k2/2
      k5=-(k2-5)
      k4=-(k2-4)
      do 10 i=1,n1
        rr=0.0
        do 101 k=1,3
101       rr=rr+d1(k,i)**2
10      dd1(i)=sqrt(rr)
      do 11 j=1,n2
        rr=0.0
        do 111 k=1,3
111        rr=rr+d2(k,j)**2
11      dd2(j)=sqrt(rr)
      do 1 i=1,n1
        do 1 j=1,n1
          dpd=dd1(i)+dd2(j)
          dmd=dd1(i)-dd2(j)
          elr=elr-cij(i,j)*(
     -     ((rc+dpd)**k5+(rc-dpd)**k5-(rc+dmd)**k5-(rc-dmd)**k5)/k5
     -     +(-dpd*(rc+dpd)**k4+dpd*(rc-dpd)**k4-dmd*(rc+dmd)**k4
     -       -dmd*(rc-dmd)**k4)/k4)/(dd1(i)*dd2(j))
1     continue
c     Densities used later will be in angstrom**3 already
      elr = elr*den*pi/(4*(k-1)*(k2-3))
      return
      end
      subroutine initgr(crm,natoms,edge)
c=====initialize the grid cavity search arrays
      dimension crm(3,natoms)
      common /cavdat/ rsph,rsphs,corlng,edgrs,edgrs2,edgrq2
      common /grid/ ncover(100,100,100),igfree(1000000),ngfree,ngrid,
     -  ng1,ng2,ng3,fltngr,fltng3,grid,e2g,e2gm,rsppe2,rspme2
      common /nmbrs/ autokc,autoa,rdtodg,boltz,pi,pi2
c     ngfree: number of uncovered grid points
c     igfree: list of uncovered gridpoints
c     ncover(i,j,k) : if > 0, the number of atoms covering that grid,
c                     if < 0, the location of that grid in igfree
      ng1=128
      ng2=ng1*ng1
      ng3=ngrid**3
      fltngr=float(ngrid)
      fltng3=float(ng3)
      grid=edge/fltngr
      e2g=-edge/2.0+edge/float(2*ngrid)
      e2gm=-edge/2.0-edge/float(2*ngrid)
      rsppe2=rsph+e2g
      rspme2=rsph-e2g
      rsphs=(rsph)**2
      do 10 i=1,ngrid
        do 10 j=1,ngrid
          do 10 k=1,ngrid
10          ncover(i,j,k)=0
c     For each atom, count the coverings of their grids
      do 11 ia=1,natoms
        imin=int((crm(1,ia)-rsppe2)/grid)+2
        imax=int((crm(1,ia)+rspme2)/grid)+1
        jmin=int((crm(2,ia)-rsppe2)/grid)+2
        jmax=int((crm(2,ia)+rspme2)/grid)+1
        kmin=int((crm(3,ia)-rsppe2)/grid)+2
        kmax=int((crm(3,ia)+rspme2)/grid)+1
c       Perform the changes
        dx=e2gm+float(imin-1)*grid-crm(1,ia)
        do 13 i=imin,imax
          dx=dx+grid
          r2x=dx*dx
          ic=i
          if (ic .lt. 1) ic=ic+ngrid
          if (ic .gt. ngrid) ic=ic-ngrid
          dy=e2gm+float(jmin-1)*grid-crm(2,ia)
          do 13 j=jmin,jmax
            dy=dy+grid
            r2y=r2x+dy*dy
            if (r2y .gt. rsphs) go to 13
            jc=j
            if (jc .lt. 1) jc=jc+ngrid
            if (jc .gt. ngrid) jc=jc-ngrid
            dz=e2gm+float(kmin-1)*grid-crm(3,ia)
            do 130 k=kmin,kmax
              dz=dz+grid
              if (r2y+dz*dz .gt. rsphs) go to 130
              kc=k
              if (kc .lt. 1) kc=kc+ngrid
              if (kc .gt. ngrid) kc=kc-ngrid
              ncover(ic,jc,kc)=ncover(ic,jc,kc)+1
130         continue
13      continue
11    continue
c     Count and collect the free sites
      ngfree=0
      do 12 i=1,ngrid
        do 12 j=1,ngrid
          do 12 k=1,ngrid
            if (ncover(i,j,k) .gt. 0) go to 12
c           Free grid (cavity) found
            indx=(i-1)*ng2+(j-1)*ng1+k
            ngfree=ngfree+1
            if (ngfree .gt. 1000000) go to 99
            ncover(i,j,k)=-ngfree
            igfree(ngfree)=indx
12    continue
      return
99    write (6,1000)
1000  format(' ***** ERROR: number of cavities in the initial ',
     -  'configuration exceeds capacity (1000000)',/,
     -  ' - recompile with increased # CV')
      stop
      end
      subroutine move(isolv,nmolec,c,crm,orient,etoto,etotn,igrid,cut2,
     -  usum,uusum,emin,emax,nemin,nemax,tesi,tex,texo,ex,nmacc,volfac)
c*****Perform a molecule displacement try
      dimension c(3,10000),crm(3,1000),orient(3,3,1000)
      real*8 tesi,etoto,etotn,usum,uusum
      dimension tesi(1000),tex(1000),texo(1000),ex(499500)
      character*80 label
      common /id/ label(2),temp,fkt,cutof,cutof2,edgea,edge,hlfedg,
     -  prsag,prsar,pressg,pressr,bparg,bparr
      common /sreal/ eming,eminr,emaxg,emaxr,vol,volg,volr,vfacg,vfacr,
     -  cut2g,cut2r,vming,vmaxg,vminr,vmaxr
      common /svint/ nmc,nmaccg,nmaccr,neming,neminr,nemaxg,nemaxr,
     -  nmcmax,nmcrep,nrecd,nmcvch,nsync,niaccg,ndaccg,nmovg,nmovr,
     -  nvming,nvmaxg,npming,npmaxg,nvminr,nvmaxr,npminr,npmaxr,
     -  nvacc,nvaccg,nvaccr,nmccnt,ncnt
      common /info/ nslv,nmolg,nmolr,natg,natr,icopvg,icopvr,nmolmx
      common /dploy/ rsi(3,10),crmi(3),orj(3,3),cedgs,rotax,delvol
      common /cavdat/ rsph,rsphs,corlng,edgrs,edgrs2,edgrq2
      common /seed/ ix
      common /option/ iop(40)
      common /nmbrs/ autokc,autoa,rdtodg,boltz,pi,pi2
      if (nmolec .gt. 0) then
        call deploy(c,crm,orient,isolv,nmolec,volfac)
c       Calculation of the total  energy
        call eism1(c,crm,tex,tesi,etotn,etoto,isolv,nmolec,cut2,volfac)
        if (etotn .le. etoto) go to 909
c       Test for acceptance
        ediff=-etotn+etoto
        scdif=exp(ediff*fkt)
        call random(1,x)
        if (scdif .gt. x) go to 909
c       Configuration rejected, restore c,crm
        is=-isolv
        call trnsfr(crm(1,isolv),crmi,3)
        if (nslv .eq. 1) go to 609
        ifat=(isolv-1)*nslv+1
        call trnsfr(c(1,ifat),rsi,3*nslv)
        go to 609
909     continue
c       Configuration accepted
        is=isolv
        call mnmxup(sngl(etotn),nmc,emin,emax,nemin,nemax)
        call getegy(texo,ex,isolv,nmolec)
        call savegy(tesi,tex,texo,ex,isolv,nmolec)
        etoto=etotn
c       Save the new rotation matrix
        if (nslv .gt. 1) call trnsfr(orient(1,1,isolv),orj,9)
        nmacc=nmacc+1
        if (igrid .gt. 1) then
          call remcov(crmi)
          call addcov(crm(1,isolv))
        end if
609     continue
        usum=usum+etoto
        uusum=uusum+etoto*etoto
      else
        is=0
        nmacc=nmacc+1
      end if
      if (mod(nmc,nmcrep) .ne. 0 .and. nmc .ne. nmcmax) return
c     Calculation of thermodynamic properties for print
      x=float(nmc+(nmc-1)/nmcvch)
      if (iop(6) .gt. 0) x=x+float(nmc-1)
      u=usum/x
      ecv=uusum/x
      acc=float(nmacc)/float(nmc)
      nemn=nemin/1000
      nemx=nemax/1000
      write (6,1015) nmc,etotn,u,ecv,emin,nemn,emax,nemx,acc,is
      return
1015  format(' N=',i7,' e=',e12.6,'  <e>=',f12.6,' <e**2>=',
     -  e12.6,' emn=',f9.2,' (',i4,'k) emx=',f9.2,' (',i4,'k)',
     -  ' a=',f5.3,' is=',i4)
      end
      subroutine indlgr(c,crm,orient,tex,texo,tesi,etoto,etotn,nmolec,
     -  cut2,indel,is,igrid,pcav,volfac)
c=====perform an insertion or deletion attempt based on grid search
      real*8 tesi,etoto,etotn
      dimension c(3,10000),crm(3,1000),orient(3,3,1000),tex(1000),
     -  texo(1000),tesi(1000)
      common /grid/ ncover(100,100,100),igfree(1000000),ngfree,ngrid,
     -  ng1,ng2,ng3,fltngr,fltng3,grid,e2g,e2gm,rsppe2,rspme2
      character*80 label
      common /id/ label(2),temp,fkt,cutof,cutof2,edgea,edge,hlfedg,
     -  prsag,prsar,pressg,pressr,bparg,bparr
      common /sreal/ eming,eminr,emaxg,emaxr,vol,volg,volr,vfacg,vfacr,
     -  cut2g,cut2r,vming,vmaxg,vminr,vmaxr
      common /svint/ nmc,nmaccg,nmaccr,neming,neminr,nemaxg,nemaxr,
     -  nmcmax,nmcrep,nrecd,nmcvch,nsync,niaccg,ndaccg,nmovg,nmovr,
     -  nvming,nvmaxg,npming,npmaxg,nvminr,nvmaxr,npminr,npmaxr,
     -  nvacc,nvaccg,nvaccr,nmccnt,ncnt
      common /info/ nslv,nmolg,nmolr,natg,natr,icopvg,icopvr,nmolmx
      common /dploy/ rsi(3,10),crmi(3),orj(3,3),cedgs,rotax,delvol
      real*8 besg,besr,etotog,etotng,etotor,etotnr
      common /energy/ besg(1000,2),texg(1000,2),texog(1000,2),
     -                exg(499500,2),
     -                besr(1000,2),texr(1000,2),texor(1000,2),
     -                exr(499500,2),
     -                etotog,etotng,etotor,etotnr
      common /option/ iop(40)
      common /nmbrs/ autokc,autoa,rdtodg,boltz,pi,pi2
      common /cavdat/ rsph,rsphs,corlng,edgrs,edgrs2,edgrq2
      dimension rno(3),crmn(3)
      call random(1,rno)
      if (indel .eq. -1) then
c       Try to *delete*
        is=rno(1)*nmolec+1
        etotn=etoto-tesi(is)
        if (igrid .gt. 1) then
          call remcov(crm(1,is))
          if (ngfree .gt. 0) then
            pcav=float(ngfree)/fltng3
          else
            pcav=1.0
          end if
        else
          pcav=1.0
        end if
        nmolec=nmolec-1
      else
c       Try to *insert*
c       Insert new molecule
        nmolec=nmolec+1
        is=nmolec
        if (nmolec .gt. 1000) then
          write (6,1000)
          call csave(10)
          stop
        end if
43      if (ngfree .eq. 0 .or. igrid .lt. 2) then
c         Do random insertion
          call ranpos(crm(1,nmolec),edge,hlfedg)
          pcav=1.0
        else
c         Find a random free grid, insert
          insx=int(rno(1)*ngfree)+1
          indx=igfree(insx)
          ix1=indx-1
          indx1=ix1/ng1
          iz=ix1-indx1*ng1
          ix=indx1/ng1
          iy=indx1-ix*ng1
          crm(1,nmolec)=e2g+ix*grid
          crm(2,nmolec)=e2g+iy*grid
          crm(3,nmolec)=e2g+iz*grid
          pcav=float(ngfree)/fltng3
        end if
c       Insert new molec
        if (nslv .gt. 1) then
          call ranort(orient(1,1,nmolec))
          do 50 k=1,3
50          crmn(k)=crm(k,nmolec)*volfac
          call putmol(nmolec,c,crmn,orient)
        end if
        tesi(nmolec)=0.d0
        call eism1(c,crm,tex,tesi,etotn,etoto,nmolec,nmolec,cut2,volfac)
      end if
      return
1000  format(' No space for the insertion')
      end
      subroutine addcov(r)
c=====update grid data when a new atom is added
      dimension r(3)
      common /cavdat/ rsph,rsphs,corlng,edgrs,edgrs2,edgrq2
      common /grid/ ncover(100,100,100),igfree(1000000),ngfree,ngrid,
     -  ng1,ng2,ng3,fltngr,fltng3,grid,e2g,e2gm,rsppe2,rspme2
      common /option/ iop(40)
c     Find the index limits for which possible change in cover can occur
      imin=int((r(1)-rsppe2)/grid)+2
      imax=int((r(1)+rspme2)/grid)+1
      jmin=int((r(2)-rsppe2)/grid)+2
      jmax=int((r(2)+rspme2)/grid)+1
      kmin=int((r(3)-rsppe2)/grid)+2
      kmax=int((r(3)+rspme2)/grid)+1
c     Perform the changes
      dx=e2gm+float(imin-1)*grid-r(1)
      do 10 i=imin,imax
        dx=dx+grid
        r2x=dx*dx
        ic=i
        if (ic .lt. 1) ic=ic+ngrid
        if (ic .gt. ngrid) ic=ic-ngrid
        dy=e2gm+float(jmin-1)*grid-r(2)
        do 10 j=jmin,jmax
          dy=dy+grid
          r2y=r2x+dy*dy
          if (r2y .gt. rsphs) go to 10
          jc=j
          if (jc .lt. 1) jc=jc+ngrid
          if (jc .gt. ngrid) jc=jc-ngrid
          dz=e2gm+float(kmin-1)*grid-r(3)
          if (kmin .gt. 0) go to 110
          kmin1=kmin+ngrid
          kmax1=ngrid
          kmin2=1
          kmax2=kmax
          go to 111
110       if (kmax .le. ngrid) go to 120
          kmin1=kmin
          kmax1=ngrid
          kmin2=1
          kmax2=kmax-ngrid
          go to 111
120       kmin2=kmin
          kmax2=kmax
          go to 112
111       do 102 kc=kmin1,kmax1
            dz=dz+grid
            if (r2y+dz*dz .gt. rsphs) go to 102
            if (ncover(ic,jc,kc) .gt. 0) go to 101
c           Grid was free
            indxx=-ncover(ic,jc,kc)
            indx=igfree(ngfree)
            igfree(indxx)=indx
            ix1=indx-1
            indx1=ix1/ng1
            iz=indx-indx1*ng1
            ix=indx1/ng1
            iy=indx1-ix*ng1+1
            ix=ix+1
            ncover(ix,iy,iz)=-indxx
            ngfree=ngfree-1
            ncover(ic,jc,kc)=1
            go to 102
101         ncover(ic,jc,kc)=ncover(ic,jc,kc)+1
102       continue
112       do 104 kc=kmin2,kmax2
            dz=dz+grid
            if (r2y+dz*dz .gt. rsphs) go to 104
            if (ncover(ic,jc,kc) .gt. 0) go to 103
c           Grid was free
            indxx=-ncover(ic,jc,kc)
            indx=igfree(ngfree)
            igfree(indxx)=indx
            ix1=indx-1
            indx1=ix1/ng1
            iz=indx-indx1*ng1
            ix=indx1/ng1
            iy=indx1-ix*ng1+1
            ix=ix+1
            ncover(ix,iy,iz)=-indxx
            ngfree=ngfree-1
            ncover(ic,jc,kc)=1
            go to 104
103         ncover(ic,jc,kc)=ncover(ic,jc,kc)+1
104       continue
10    continue
      return
      end
      subroutine remcov(r)
c=====update grid data when an atom is removed
      dimension r(3)
      common /cavdat/ rsph,rsphs,corlng,edgrs,edgrs2,edgrq2
      common /grid/ ncover(100,100,100),igfree(1000000),ngfree,ngrid,
     -  ng1,ng2,ng3,fltngr,fltng3,grid,e2g,e2gm,rsppe2,rspme2
      common /option/ iop(40)
c     Find the index limits for which possible change in cover can occur
      imin=int((r(1)-rsppe2)/grid)+2
      imax=int((r(1)+rspme2)/grid)+1
      jmin=int((r(2)-rsppe2)/grid)+2
      jmax=int((r(2)+rspme2)/grid)+1
      kmin=int((r(3)-rsppe2)/grid)+2
      kmax=int((r(3)+rspme2)/grid)+1
c     Perform the changes
      dx=e2gm+float(imin-1)*grid-r(1)
      do 10 i=imin,imax
        dx=dx+grid
        r2x=dx*dx
        ic=i
        if (ic .lt. 1) ic=ic+ngrid
        if (ic .gt. ngrid) ic=ic-ngrid
        dy=e2gm+float(jmin-1)*grid-r(2)
        do 10 j=jmin,jmax
          dy=dy+grid
          r2y=r2x+dy*dy
          if (r2y .gt. rsphs) go to 10
          jc=j
          if (jc .lt. 1) jc=jc+ngrid
          if (jc .gt. ngrid) jc=jc-ngrid
          dz=e2gm+float(kmin-1)*grid-r(3)
          if (kmin .gt. 0) go to 110
          kmin1=kmin+ngrid
          kmax1=ngrid
          kmin2=1
          kmax2=kmax
          go to 111
110       if (kmax .le. ngrid) go to 120
          kmin1=kmin
          kmax1=ngrid
          kmin2=1
          kmax2=kmax-ngrid
          go to 111
120       kmin2=kmin
          kmax2=kmax
          go to 112
111       do 102 kc=kmin1,kmax1
            dz=dz+grid
            if (r2y+dz*dz .gt. rsphs) go to 102
            if (ncover(ic,jc,kc) .gt. 1) go to 101
c           Grid became free
            ngfree=ngfree+1
            igfree(ngfree)=(ic-1)*ng2+(jc-1)*ng1+kc
            ncover(ic,jc,kc)=-ngfree
            go to 102
101         ncover(ic,jc,kc)=ncover(ic,jc,kc)-1
102       continue
112       do 104 kc=kmin2,kmax2
            dz=dz+grid
            if (r2y+dz*dz .gt. rsphs) go to 104
            if (ncover(ic,jc,kc) .gt. 1) go to 103
c           Grid became free
            ngfree=ngfree+1
            igfree(ngfree)=(ic-1)*ng2+(jc-1)*ng1+kc
            ncover(ic,jc,kc)=-ngfree
            go to 104
103         ncover(ic,jc,kc)=ncover(ic,jc,kc)-1
104       continue
10    continue
      if (ngfree .le. 1000000) return
      write (6,1000) ngfree
1000  format(' Number of free gridpoints=',i8,' exceeds the capacity ',
     -  'of the program (1000000)',/,' - recompile with increased # CV')
      stop
      end
      subroutine crstr(iunit)
c*****Save all the common blocks
      character*80 c2
      real*8 d2,d10,d12,d22
      common /seed/ n1
      common /id/ c2(40),x1(13)
      common /accu/ d2(15)
      common /sreal/ x18(4),x3(11)
      common /svint/ n3(29)
      common /info/ n4(8)
      common /grid/ n5(100,100,100),n6(1000000),n7(5),x4(7)
      common /coord/ x5(7,10),n8(10)
      common /dploy/ x6(3,10),x7(15)
      common /ctotal/ x8(12,10000),x9(24,1000)
      common /energy/ d10(4,1000),x11(8,1000),d12(4),x13(4,499500)
      common /cavdat/ x14(6)
      common /cavprb/ n9(2,1000),x17(4)
      common /ecell/ x15(81)
      common /option/ n10(40)
      common /spcpar/ x16(9)
      common /genslv/ x19(4,10),x20(3,10,10)
      common /cntrl/ x21(4,100)
      common /distr/ d22(2,210),n11(22,210),x24(20,210),x25(31),n12(6)
      rewind iunit
      read (iunit)
     -   n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,
     -   x1,d2,x3,x4,x5,x6,x7,x8,x9,d10,x11,d12,x13,x14,x15,x16,x17,x18,
     -   x19,x20,x21,d22,x23,x24,x25
      return
      end
      subroutine csave(iunit)
c*****Restore all the common blocks
      character*80 c2
      real*8 d2,d10,d12,d22
      common /seed/ n1
      common /id/ c2(40),x1(13)
      common /accu/ d2(15)
      common /sreal/ x18(4),x3(11)
      common /svint/ n3(29)
      common /info/ n4(8)
      common /grid/ n5(100,100,100),n6(1000000),n7(5),x4(7)
      common /coord/ x5(7,10),n8(10)
      common /dploy/ x6(3,10),x7(15)
      common /ctotal/ x8(12,10000),x9(24,1000)
      common /energy/ d10(4,1000),x11(8,1000),d12(4),x13(4,499500)
      common /cavdat/ x14(6)
      common /cavprb/ n9(2,1000),x17(4)
      common /ecell/ x15(81)
      common /option/ n10(40)
      common /genslv/ x19(4,10),x20(3,10,10)
      common /spcpar/ x16(9)
      common /cntrl/ x21(4,100)
      common /distr/ d22(2,210),n11(22,210),x24(20,210),x25(31),n12(6)
      rewind iunit
      write (iunit)
     -   n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,
     -   x1,d2,x3,x4,x5,x6,x7,x8,x9,d10,x11,d12,x13,x14,x15,x16,x17,x18,
     -   x19,x20,x21,d22,x23,x24,x25
      return
      end
      subroutine echopr
c*****Print all the parameters of the run
      common /option/ iop(40)
      character*80 label
      common /id/ label(2),temp,fkt,cutof,cutof2,edgea,edge,hlfedg,
     -  prsag,prsar,pressg,pressr,bparg,bparr
      common /spcpar/qqspc,qqspc2,qqspc4,alj,blj,alj6,alj66,blj12,blj122
      common /genslv/ c6(10),c12(10),cc6(10,10),cc12(10,10),qqg(10,10),
     -                siglj(10),epslj(10)
      common /sreal/ eming,eminr,emaxg,emaxr,vol,volg,volr,vfacg,vfacr,
     -  cut2g,cut2r,vming,vmaxg,vminr,vmaxr
      common /svint/ nmc,nmaccg,nmaccr,neming,neminr,nemaxg,nemaxr,
     -  nmcmax,nmcrep,nrecd,nmcvch,nsync,niaccg,ndaccg,nmovg,nmovr,
     -  nvming,nvmaxg,npming,npmaxg,nvminr,nvmaxr,npminr,npmaxr,
     -  nvacc,nvaccg,nvaccr,nmccnt,ncnt
      common /dploy/ rsi(3,10),crmi(3),orj(3,3),cedgs,rotax,delvol
      common /info/ nslv,nmolg,nmolr,natg,natr,icopvg,icopvr,nmolmx
      common /grid/ ncover(100,100,100),igfree(1000000),ngfree,ngrid,
     -  ng1,ng2,ng3,fltngr,fltng3,grid,e2g,e2gm,rsppe2,rspme2
      common /coord/ relcrd(3,10),cslv(3,10),qslv(10),islv(10)
      common /cavdat/ rsph,rsphs,corlng,edgrs,edgrs2,edgrq2
      real*8 angcrg,angcrr
      common /distr/ angcrg(210),angcrr(210),nrg(210),nrr(210),
     -  rdbg(10,210),rdbr(10,210),angcbg(10,210),angcbr(10,210),
     -  ncnfbg(10),ncnfbr(10),rgrid,rgridi,rgrdig,rgrdir,rdmax2,rdmx2g,
     -  rdmx2r,rdm2ng,rdm2nr,rgring,rgrinr,nrgrid,icolcr,nconfg,nconfr,
     -  lenblk,nblk
      character*4 iel
      common /iname/ iel(100)
      common /nmbrs/ autokc,autoa,rdtodg,boltz,pi,pi2
      if (iop(1) .eq. 1) then
        write (6,011) nmolg,nmolr
      else
        write (6,010)
        if (iop(4) .eq. 1) write (6,040)
        if (iop(4) .eq. 2) write (6,012)
      end if
      if (iop(3) .eq. 0) write (6,100)
      if (iop(3) .eq. 1) write (6,101) prsag,prsar
      if (iop(3) .eq. 2) write (6,102) bparg,bparr
      if (iop(11) .eq. 0) write (6,113)
      if (iop(11) .eq. 1) write (6,114)
      if (iop(11) .eq. 2) write (6,126) alj,blj
      if (iop(11) .eq. 3) write (6,127) alj,blj
      if (iop(11) .eq. 4) write (6,111)
      if (iop(11) .eq. 5) write (6,116)
      if (iop(11) .eq. 6) write (6,112)
      if (iop(6) .eq. 0) write (6,060)
      if (iop(6) .eq. 1) write (6,061)
      if (iop(6) .eq. 2) write (6,062) ngrid,rsph
      if (iop(17) .gt. 0) then
        rmax=nrgrid*rgrid
        write (6,171) rgrid,rmax
        if (iop(17) .gt. 1) write (6,172) icolcr
        if (iop(17) .gt. 2) write (6,174)
        write (6,173) lenblk
      end if
      write (6,103)
      qsum=0.0
      do 200 i=1,nslv
        qsum=qsum+qslv(i)
        iam=islv(i)
        cx=cslv(1,i)
        cy=cslv(2,i)
        cz=cslv(3,i)
c        if (iop(11) .le. 4 .or. c6(i)*c12(i) .eq. 0.0) then
        if (iop(11) .le. 4 .or. siglj(i)*epslj(i) .eq. 0.0) then
          write (6,106) i,iel(iam),cx,cy,cz,qslv(i),c6(i),c12(i)
        else
          sig=(c12(i)/c6(i))**(1.0/6.0)
          eps=c6(i)**2/(4.0*c12(i))
          write (6,104) i,iel(iam),cx,cy,cz,qslv(i),c6(i),c12(i),eps,sig
        end if
200   continue
      write (6,105) qsum
      rota=rotax*rdtodg
      write (6,1013) temp,edgea,cutof
      if (iop(8) .ne. 0) write (6,117) corlng
      write (6,107) cedgs,rota,delvol
      write (6,108) nmcmax,nmcrep,nrecd,nmcvch,nsync
c      if (1000000 .ge. 2**15-1) then
c        write (6,130)
c        stop
c      end if
      if (nslv .gt. 10 .or. ngrid .gt. 100 .or.
     -  nmolg .gt. 1000 .or. nmolr .gt. 1000 .or.
     -  nmolg*(nmolg-1)/2 .gt. 499500 .or. nmolg*nslv .gt. 10000 .or.
     -  nmolr*(nmolr-1)/2 .gt. 499500 .or. nmolr*nslv .gt. 10000) then
        write (6,131)
        stop
      end if
999   return
010   format(' Input configurations read from unit 13')
011   format(' New initial configuration was created with',i5,' and ',
     -  i5,' molecules')
012   format(' Number of molecules adjusted from input on unit 5')
040   format(' Alphanumeric input from unit 13')
060   format(' No particle exchange is done')
061   format(' Insertion algorithm: random in both systems')
062   format(' Insertion algorithm: cavity-biased with grid search in',
     -  ' g system, random in the r system',/,
     -  5x,' Number of grids along each axis=',i5,/,
     -  5x,' Cavity radius=',f7.3,' A')
100   format(' Gibbs ensemble simulation')
101   format(' Two parallel (T,P,N) ensemble simulations at pressures=',
     -  2f10.5,' atm')
102   format(' Two parallel (T,V,mu) ensemble simulations at b param=',
     -  2f10.5)
103   format(/,' Nuclear geometry:',//,' no',17x,'x',16x,'y',16x,'z',
     -  16x,'q',/)
104   format(i3,2x,a4,2x,4(2x,f15.8),2f12.2,' eps=',f6.4,' sig=',f6.3)
105   format(' The charge sum=',f10.5)
106   format(i3,2x,a4,2x,4(2x,f15.8),2f12.2)
107   format(' Shift interval   = ',f10.4,' A',/,' Rotation interval=',
     -  f10.4,' deg ',/,' Volume change interval=',f10.4,' A**3')
108   format(/,' Number of MC steps to run=',i9,/,
     -  ' Result interval=',i7,' MC steps',/,
     -  ' Record interval=',i7,' MC steps',/,
     -  ' Volume change attempt frequency=',i7,' MC steps',/,
     -  ' Frequency of c.o.m. synchronization=',i9,' MC steps')
1013  format(' Absolute temperature=',f10.4,' Kelvin',/,
     -  ' Size of the original unit cells=',f10.5,' A',/,
     -  ' Potential cutoff=',f12.5,' A')
111   format(' Lennard-Jones fluid')
112   format(' General 12-6-1 fluid')
113   format(' MCY water')
114   format(' ST2 water')
116   format(' General 12-6 fluid')
117   format(' Long range correction*vol/n**2=',e15.5)
126   format(' SPC/TIPS water  alj=',f15.6,'  blj=',f15.6)
127   format(' TIPS2/TIPS4P water  alj=',f15.6,'  blj=',f15.6)
130   format(' Too many free grid points allowed - recompile with',
     -  ' the array ngrid as integer*4')
131   format(' Input data for nmolec, nslv, or ngrid exceeds program ',
     -  'dimensions - recompile with increased # MO, # SV or # GX')
171   format(' Radial distribution function will be calculated with',
     -  ' a gridsize of ',f6.2,' A to a distance of ',f8.2,' A')
172   format(' Orientational correlation based on the ',i2,'-th axis',
     -  ' will also be calculated')
173   format(' Radial distribution function control function block',
     -  ' size=',i8,' Monte Carlo steps') 
174   format(' The absolute value of the orientational correlations',
     -  ' will be accumulated')
      end
      subroutine eism0(c,crm,ex,tesi,esm,nmolec,cut2,volfac,
     -  irdf,nr,angcr,or,nconf,icolcr,rgridi,rdmax2)
c*****Compute the total energy esm of a configuration from scratch
      real*8 esm,tesi,angcr
      dimension tesi(1000),tex(1000),ex(499500),c(3,10000),crm(3,1000),
     -  nr(210),angcr(210),or(3,3,1000)
      common /texz/ tex0(1000)
      common /nmbrs/ autokc,autoa,rdtodg,boltz,pi,pi2
      dimension rrj(1000),ori(3)
      esm=0.d0
      do i=1,nmolec
        tesi(i)=0.d0
      end do
      if (nmolec .eq. 1) return
      do i=2,nmolec
        i1=i-1
        do j=1,i1
          tex(j)=engclc(c,crm,i,j,ipbc,cut2,volfac,rr)
          esm=esm+dble(tex(j))
          rrj(j)=rr
        end do
        call savegy(tesi,tex,tex0,ex,i,i)
        if (irdf .gt. 0) then
          do k=1,3
            ori(k)=or(k,icolcr,i)
          end do
          do j=1,i1
            if (rrj(j) .le. rdmax2) then
              ir=sqrt(rrj(j))*rgridi
              nr(ir)=nr(ir)+1
              if (irdf .gt. 1) then
c               Accumulate orientational correlations
                cr=0.0
                do k=1,3
                  cr=cr+or(k,icolcr,j)*ori(k)
                end do
                if (irdf .gt. 2) cr=abs(cr)
                angcr(ir)=angcr(ir)+cr
              end if
            end if
          end do
        end if
      end do
      if (irdf .gt. 0) nconf=nconf+1
      return
      end
      subroutine eism1(c,crm,tex,tesi,esmn,esmo,is,nmolec,cut2,volfac)
c*****Update the total energy after having moved the is-th molecule
c*****old value: esmo, new value:esmn.
      real*8 esmn,esmo,esmnew,tesi
      dimension tex(1000),tesi(1000),c(3,10000),crm(3,1000)
      common /nmbrs/ autokc,autoa,rdtodg,boltz,pi,pi2
      esmnew=0.d0
      tex(is)=0.0
      i1=is-1
      i2=is+1
      if (i1 .ne. 0) then
        do i=1,i1
          tex(i)=engclc(c,crm,is,i,ipbc,cut2,volfac,rr)
          esmnew=esmnew+dble(tex(i))
        end do
      end if
      if (i2 .le. nmolec) then
        do i=i2,nmolec
          tex(i)=engclc(c,crm,i,is,ipbc,cut2,volfac,rr)
          esmnew=esmnew+dble(tex(i))
        end do
      end if
      esmn=esmo+(esmnew-tesi(is))
      return
      end
      subroutine deploy(c,crm,orient,isolv,nmolec,volfac)
c*****Translate and rotate the isolv-th molecule
      dimension c(3,10000),crm(3,1000),orient(3,3,1000)
      character*80 label
      common /id/ label(2),temp,fkt,cutof,cutof2,edgea,edge,hlfedg,
     -  prsag,prsar,pressg,pressr,bparg,bparr
      common /coord/ relcrd(3,10),cslv(3,10),qslv(10),islv(10)
      common /nmbrs/ autokc,autoa,rdtodg,boltz,pi,pi2
      common /info/ nslv,nmolg,nmolr,natg,natr,icopvg,icopvr,nmolmx
      common /dploy/ rsi(3,10),crmi(3),orj(3,3),cedgs,rotax,delvol
      dimension rno(5),rsj(3,10),ch(3),crmfac(3)
      call random(5,rno)
c     Save the old coordinates, center of masses
      call trnsfr(crmi,crm(1,isolv),3)
      if (nslv .gt. 1) then
        ifat=(isolv-1)*nslv+1
        call trnsfr(rsi,c(1,ifat),3*nslv)
c       Generate rotation
        iaxis=rno(5)*3+1
        phi=(rno(4)-0.5)*rotax
        do k=1,3
          crmfac(k)=crmi(k)*volfac
        end do
        call rotate(rsj,orient(1,1,isolv),orj,crmfac,phi,nslv,iaxis)
      end if
c     Generate translation
c     Uniform cubic move
      do k=1,3
        ch(k)=(rno(k)-0.5)*cedgs
      end do
      do k=1,3
        crm(k,isolv)=crm(k,isolv)+ch(k)
      end do
      call cellfn(rno(1),rno(2),rno(3),crm(1,isolv),crm(2,isolv),
     -  crm(3,isolv),hlfedg,ipbc)
      if (ipbc .gt. 1) then
        do k=1,3
          crm(k,isolv)=crm(k,isolv)+rno(k)
          ch(k)=ch(k)+rno(k)
        end do
      end if
      if (nslv .eq. 1) return
      do i=1,nslv
        do k=1,3
          c(k,i+ifat-1)=rsj(k,i)+ch(k)*volfac
        end do
      end do
      return
      end
      subroutine volcha(vol,voln,cold,cnew,crm,nmol,nat,volfo,volfn)
c*****Calculates new atomic coordinates for a volume change
      dimension cold(3,10000),cnew(3,10000),crm(3,1000)
c     volfo, volfn are the scaling factors corresponding to the old and new vols
      dimension delcrm(3,1000)
c     Loops are written with vectorizability in mind.
      volfn=(voln/vol)**(1.0/3.0)
c     delcrm will contain the molecular shift vectors due to the volume change
      if (nat .gt. 1) then
        do im=1,nmol
          delcrm(1,im)=crm(1,im)*(volfn-volfo)
          delcrm(2,im)=crm(2,im)*(volfn-volfo)
          delcrm(3,im)=crm(3,im)*(volfn-volfo)
        end do
        do ia=1,nat
          do im=1,nmol
            cnew(1,(im-1)*nat+ia)=cold(1,(im-1)*nat+ia)+delcrm(1,im)
            cnew(2,(im-1)*nat+ia)=cold(2,(im-1)*nat+ia)+delcrm(2,im)
            cnew(3,(im-1)*nat+ia)=cold(3,(im-1)*nat+ia)+delcrm(3,im)
          end do
        end do
      end if
      return
      end
      subroutine rotate(csb,ori,orj,rmass,phi,natm,iaxis)
c*****Rotate a molecule by phi around iaxis
      common /coord/ relcrd(3,10),cslv(3,10),qslv(10),islv(10)
      dimension rmass(3),rot(3,3),ori(3,3),orj(3,3),csb(3,10)
      csp=cos(phi)
      snp=sin(phi)
      i1=mod(iaxis,3)+1
      i2=mod(iaxis+1,3)+1
      i3=mod(iaxis+2,3)+1
c     rot: rotation matrix of the change
      do 15 k=1,3
        do 15 l=1,3
15        rot(k,l)=0.0
      rot(i1,i1)=csp
      rot(i1,i2)=snp
      rot(i2,i1)=-snp
      rot(i2,i2)=csp
      rot(i3,i3)=1.0
c     orj: rotation matrix transforming the local water to the actual
      do 16 i=1,3
        do 16 j=1,3
        rr=0.0
        do 17 k=1,3
17        rr=rr+rot(i,k)*ori(k,j)
16      orj(i,j)=rr
c     Make sure that the new rotation matrix is orthogonal
      call vprod(orj,1,2,3)
      call vprod(orj,3,1,2)
      call mnorm(orj)
c     Perform the transformation of the coordinates and the lone p.
      do 10 i=1,natm
        do 20 j=1,3
          rr=0.0
          do 21 k=1,3
21          rr=rr+orj(j,k)*relcrd(k,i)
20        csb(j,i)=rr+rmass(j)
10    continue
      return
      end
      subroutine ranpos(crm,edge,hlfedg)
c*****Create a random location in the simulation cell
c     crm is in reference volume scale
      dimension crm(3),ran(3)
      call random(3,ran)
c     Simple cubic
      crm(1)=-hlfedg+edge*ran(1)
      crm(2)=-hlfedg+edge*ran(2)
      crm(3)=-hlfedg+edge*ran(3)
      return
      end
      subroutine putmol(nmol,c,crm,orient)
c*****Create a molecule
c     crm is in current volume scale
      dimension c(3,10000),crm(3),orient(3,3,1000)
      common /info/ nslv,nmolg,nmolr,natg,natr,icopvg,icopvr,nmolmx
      common /coord/ relcrd(3,10),cslv(3,10),qslv(10),islv(10)
      do 20 i=1,nslv
        do 20 j=1,3
          rr=0.0
          do 21 k=1,3
21          rr=rr+orient(j,k,nmol)*relcrd(k,i)
20        c(j,(nmol-1)*nslv+i)=crm(j)+rr
      return
      end
      subroutine genmol(nmol1,nmol2,nmol,c,crm,or,edge,hlfedg,volfac)
c*****Create new molecules from nmol1+1 to nmol2
      dimension c(3,10000),crm(3,1000),or(3,3,1000),crmi(3)
      if (nmol1 .eq. 0) write(6,1001) nmol2
      if (nmol1 .gt. 0) write(6,1006) nmol1,nmol2
      im1=nmol1+1
      do 2624 im=im1,nmol2
        call ranpos(crm(1,im),edge,hlfedg)
        call ranort(or(1,1,im))
        do 2600 k=1,3
2600      crmi(k)=crm(k,im)*volfac
2624    call putmol(im,c,crmi,or)
      nmol=nmol2
      return
1001  format(1x,i4,' molecules created')
1006  format(' Number of molecules=',i4,' extended to ',i4)
      end
      subroutine ranort(r)
      dimension r(3,3),rn(3)
c*****This subroutine prepares a random rotation matrix
c     from the 3 Euler angles
c     Eq(4-47) of Goldstein (r=a)
      common /nmbrs/ autokc,autoa,rdtodg,boltz,pi,pi2
      call random(3,rn)
      fir=rn(1)*pi2
      psr=rn(3)*pi2
      cth=-1.0+2.0*rn(2)
      sth=sqrt(1.0-cth*cth)
      sfi=sin(fir)
      cfi=cos(fir)
      sps=sin(psr)
      cps=cos(psr)
      r(1,1)=cps*cfi-cth*sfi*sps
      r(2,1)=-sps*cfi-cth*sfi*cps
      r(3,1)=sth*sfi
      r(1,2)=cps*sfi+cth*cfi*sps
      r(2,2)=-sps*sfi+cth*cfi*cps
      r(3,2)=-sth*cfi
      r(1,3)=sth*sps
      r(2,3)=sth*cps
      r(3,3)=cth
      return
      end
      subroutine grelcd(cslv,relcrd,islv,nslv)
c*****Obtain the com-centered coordinates of the molecule
      dimension cslv(3,nslv),relcrd(3,nslv),islv(nslv),cmcrd(3)
      call cofms(cslv,cmcrd,islv,1,nslv)
      do 2500 i=1,nslv
        do 2500 k=1,3
2500      relcrd(k,i)=cslv(k,i)-cmcrd(k)
      return
      end
      subroutine findor(c,crm,orient,nmolec,nslv,islv,relcrd,volfac,
     -  e2,it)
c*****Obtain the rotation matrices for each molecule from c
      character*8 it
      dimension  c(3,10000),crm(3,1000),orient(3,3,1000),relcrd(3,nslv),
     -    islv(nslv),crmi(3),idel(1000)
      if (nmolec .eq. 0) return
c     Determine the rotation matrix for the solvent molecules
      ndel=0
      do 2601 im=1,nmolec
        call getort(c(1,(im-1)*nslv+1),crmi,islv,relcrd,
     -      nslv,orient(1,1,im))
        d=0.0
        do 2602 k=1,3
          crmi(k)=crmi(k)/volfac
2602      d=d+abs(crm(k,im)-crmi(k))
        if (d .gt. 0.001) write (6,1001) im,d,crmi,(crm(k,im),k=1,3),it
        do 2603 k=1,3
2603      crm(k,im)=crmi(k)
        call cellfn(cc1,cc2,cc3,crm(1,im),crm(2,im),crm(3,im),e2,ip)
        if (ip .gt. 1) then
          ndel=ndel+1
          idel(ndel)=im
        end if
2601  continue
      if (ndel .gt. 0) then
        do 2620 j=1,ndel
          i=ndel-j+1
          write (6,1000) idel(i)
          call trnsfr(crm(1,idel(i)),crm(1,nmolec),3)
          call trnsfr(c(1,(idel(i)-1)*nslv+1),
     -                c(1,(nmolec-1)*nslv+1),3*nslv)
          call trnsfr(orient(1,1,idel(i)),orient(1,1,nmolec),9)
2620      nmolec=nmolec-1
      end if
      return
1000  format(' Molecule ',i3,' was out of the box - removed')
1001  format(' ----- COM error for molecule',i4,' =',f10.5,' new=',
     -  3f10.5,' old=',3f10.5,' (',a8,')')
      end
      subroutine getort(c,crm,icl,rlc,n,orient)
c#    MMC routine 018 lstmod: before 07/15/85
c*****Obtain orientation of a molecule w.r.t. its local representaton
      dimension c(3,n),crm(3),icl(n),rlc(3,n),orient(3,3)
      dimension a(3,3),b(3,3)
      call cofms(c,crm,icl,1,n)
      call trnsfr(a,rlc,9)
      do 10 i=1,3
        do 10 k=1,3
10        b(k,i)=c(k,i)-crm(k)
      call ormat(orient,a,b,n)
      return
      end
      subroutine ormat(orm,a,b,natoms)
c#    MMC routine 022 lstmod: 12/08/86
c     3-atom part originally written by P.K. Mehrotra.
c*****Compute the rotation matrix orm  by solving the following equation
c      b=orm.a
c     ++++++++++++++++++++++++++++++++++++++
c     a - coordinates of the molecule before rotation ( local system)
c     b - coordinates of the molecule after rotation (global syetm).
c         important - only rotation , no translation.
c     ++++++++++++++++++++++++++++++
c     It is assumed that the coordinates of the molecule are given
c     columnwise, i.e. the coordinates of the atom 1 constitute
c     the first column, coordinates of the atom 2 constitute the
c     the second column, etc.
c     The elements of a and b are preserved.
      dimension a(3,3),b(3,3),orm(3,3),d(3,3),e(3,3)
c     choose the first atom as the origin
      if (natoms .eq. 1) go to 40
c     Generate three specific orthogonal vectors in lab frame
      if (natoms .eq. 2) go to 20
      do 10 i=1,3
        do 10 k=1,3
          d(k,i)=a(k,i)-a(k,1)
10        e(k,i)=b(k,i)-b(k,1)
      call vprod(d,2,3,1)
      call vprod(d,1,2,3)
      call vprod(e,2,3,1)
      call vprod(e,1,2,3)
c     Check for colinearity of the atoms
      dsum=0.0
      do 11 k=1,3
11      dsum=dsum+abs(d(k,1))
      if (dsum .gt. 1.e-6) go to 30
c     Switch to two-atom algorithm
      write (6,1000)
c     Diatomic
20    dsum=0.0
      do 21 k=1,3
        d(k,1)=a(k,2)-a(k,1)
        dsum=dsum+abs(d(k,1))
        e(k,1)=b(k,2)-b(k,1)
        d(k,2)=0.0
21      e(k,2)=0.0
c     Create the second vector as perpendicular to the bond
c     if the two atoms coincide,generate unit matrix
      if (dsum .gt. 1.e-6) go to 25
      write (6,1002)
      go to 40
25    do 22 k=1,3
        if (abs(d(k,1)) .gt. 1.e-7) go to 22
c       Zero component found
        d(k,2)=1.0
        go to 23
22    continue
c     No zero component
      d(2,2)=1.0
      d(3,2)=-d(2,1)/d(3,1)
23    call vprod(d,1,2,3)
      do 24 k=1,3
        if (abs(e(k,1)) .gt. 1.e-7) go to 24
c       Zero component found
        e(k,2)=1.0
        go to 26
24    continue
c     No zero component
      e(2,2)=1.0
      e(3,2)=-e(2,1)/e(3,1)
26    call vprod(e,1,2,3)
30    call mnorm(d)
      call mnorm(e)
c     Now, d=orm*e, thus orm=d*inv(e) and
c     the inverse of an orthonormal matrix is its transpose
      do 35 i=1,3
        do 35 j=1,3
          sum=0.0
          do 36 k=1,3
36        sum=sum+e(i,k)*d(j,k)
35    orm(i,j)=sum
      return
40    do 41 i=1,3
        do 42 j=1,3
42        orm(i,j)=0.0
41      orm(i,i)=1.0
      return
1000  format(' ----- WARNING: the first three solute atoms are on',
     -  ' the same line - two-atom algorithm will be used',/,
     -  22x,'check if solute atoms should be rearranged')
1002  format(' ----- WARNING: the two atoms of a diatomic solute ',
     -  'coincide - unit matrix will be used as orientation matrix')
      end
      subroutine cofms(csa,rmass,ian,ifat,natm)
c*****Compute center of mass of a molecule
      dimension csa(3,10),rmass(3),ian(10)
      common /atwght/ avgd,aw(100)
      fnumx=0.
      fnumy=0.
      fnumz=0.
      den=0.
      ij=ifat -1
      do i=1,natm
        ij=ij+1
        iam=ian(ij)
        den=den+aw(iam)
        fnumx=fnumx+aw(iam)*csa(1,i)
        fnumy=fnumy+aw(iam)*csa(2,i)
        fnumz=fnumz+aw(iam)*csa(3,i)
      end do
      rmass(1)=fnumx/den
      rmass(2)=fnumy/den
      rmass(3)=fnumz/den
      return
      end
      subroutine random(n,rno)
c*****Congruential random number generator, Forsythe's constants
      common /seed/ ix
      dimension rno(n)
      do 50 i=1,n
        iy=ix*314159269+453806245
c       Eliminate bits over 31
        iy=ibclr(iy,31)
        rno(i)=float(iy)/2.1474838e+09
c       rno(i)=float(iy)/2.1474836e+09
        ix=iy
50    continue
      return
      end
      subroutine ewwint
c*****This subroutine initializes the necessary variables and
c     computes certain expressions  for the water potentials
      common /cnew/ qq,dqq,a1,a2,a3,a4,b1,b2,b3,b4,ab1,ab2,ab3,ab4
      common /ch1dat/ p,q2,dq2,qq2,c1,c2,c3,c4,c5,e1,e2,e3,e4,e5,
     - ce1,ce2,ce3,ce4,ce5
      common /st2par/ qww(5,5),rl,rl2,ru2,r3ul,rul3i,sg6,eps4s
      common /spcpar/qqspc,qqspc2,qqspc4,alj,blj,alj6,alj66,blj12,blj122
      common /genslv/ c6(10),c12(10),cc6(10,10),cc12(10,10),qqg(10,10),
     -                siglj(10),epslj(10)
      common /info/ nslv,nmolg,nmolr,natg,natr,icopvg,icopvr,nmolmx
      common /coord/ relcrd(3,10),cslv(3,10),qslv(10),islv(10)
      common /nmbrs/ autokc,autoa,rdtodg,boltz,pi,pi2
      common /option/ iop(40)
      equivalence (iop11,iop(11))
      go to (110,120,120,140,150,150),iop11
c     Initialize combined ewwmcy parameters
      qq=qq*autokc*autoa
      dqq=dqq*autokc*autoa
      ab1=a1*b1
      ab2=a2*b2
      ab3=a3*b3
      ab4=a4*b4
      return
c     ST2: the water atoms are in the order o,h1,h2,lp1,lp2
110   rl=2.0160
      ru=3.1287
      rl2=rl*rl
      ru2=ru*ru
      r3ul=3.0*ru - rl
      rul3i=1.0/(ru - rl)**3
      sigma=3.10
      eps=0.07575
      sg6=sigma**6
      eps4s=4.0*eps*sg6
      do i=1,5
        do j=1,5
          qww(i,j)=qslv(i)*qslv(j)*autokc*autoa
        end do
      end do
      return
c     Initialize SPC/TIPS ...   constants
120   alj6=alj
      blj12=blj
      alj66=6.0*alj6
      blj122=12.0*blj12
      qqspc=qslv(2)**2*autokc*autoa
      qqspc2=2.0*qqspc
      qqspc4=4.0*qqspc
      return
140   continue
c     LJ fluid
      return
c     Initialize general solvent
150   do i=1,nslv
        do j=1,nslv
          sss=(siglj(i)+siglj(j))/2.
          eee=sqrt(epslj(i)*epslj(j))
          cc6(i,j)=4.*eee*sss**6
          cc12(i,j)=4.*eee*sss**12
c          cc6(i,j)=sqrt(c6(i)*c6(j))
c          cc12(i,j)=sqrt(c12(i)*c12(j))
          qqg(i,j)=qslv(i)*qslv(j)*autokc*autoa
        end do
      end do
      return
      end
      function engclc(c,crm,i,j,ipbc,cut2,volfac,rr)
c*****Extract the atomic coordinates and call the appropriate eng fct.
      dimension c(3,10000),crm(3,1000)
      character*80 label
      common /id/ label(2),temp,fkt,cutof,cutof2,edgea,edge,hlfedg,
     -  prsag,prsar,pressg,pressr,bparg,bparr
      common /info/ nslv,nmolg,nmolr,natg,natr,icopvg,icopvr,nmolmx
      common /option/ iop(40)
      dimension z(3),cc(3),b(3,10)
      equivalence (iop11,iop(11))
      do 100 k=1,3
100     z(k)=crm(k,j)-crm(k,i)
      call cellfn(cc(1),cc(2),cc(3),z(1),z(2),z(3),hlfedg,ipbc)
      if (ipbc .gt. 1) then
        do 110 k=1,3
110       z(k)=z(k)+cc(k)
      end if
      rr=z(1)*z(1)+z(2)*z(2)+z(3)*z(3)
      if (rr .gt. cut2) then
        engclc=0.0
        return
      end if
c     Skip evaluation if r56 is lt 1.5 angstroms (2.835 a.u.)
c     if (rr .le. 2.8350) then
c       engclc=1.d0+0.01d0*exp(2.8350-rr)
c       return
c     end if
      j1=(j-1)*nslv+1
      call trnsfr(b,c(1,j1),3*nslv)
      if (ipbc .gt. 1) then
        do 120 is=1,nslv
          do 120 k=1,3
120         b(k,is)=b(k,is)+cc(k)*volfac
      end if
      i1=(i-1)*nslv+1
      go to (11,12,13,14,15,16),iop11
      engclc=ewwmcy(
     -  c(1,i1),c(2,i1),c(3,i1),c(1,i1+1),c(2,i1+1),c(3,i1+1),
     -  c(1,i1+2),c(2,i1+2),c(3,i1+2),c(1,i1+3),c(2,i1+3),c(3,i1+3),
     -  b(1,1),b(2,1),b(3,1),b(1,2),b(2,2),b(3,2),
     -  b(1,3),b(2,3),b(3,3),b(1,4),b(2,4),b(3,4))
      return
11    engclc=ewwst2(c(1,i1),b)
      return
12    engclc=ewwspc(
     -  c(1,i1),c(2,i1),c(3,i1),c(1,i1+1),c(2,i1+1),c(3,i1+1),
     -  c(1,i1+2),c(2,i1+2),c(3,i1+2),
     -  b(1,1),b(2,1),b(3,1),b(1,2),b(2,2),b(3,2),
     -  b(1,3),b(2,3),b(3,3))
      return
13    engclc=ewwtip(
     -  c(1,i1),c(2,i1),c(3,i1),c(1,i1+1),c(2,i1+1),c(3,i1+1),
     -  c(1,i1+2),c(2,i1+2),c(3,i1+2),c(1,i1+3),c(2,i1+3),c(3,i1+3),
     -  b(1,1),b(2,1),b(3,1),b(1,2),b(2,2),b(3,2),
     -  b(1,3),b(2,3),b(3,3),b(1,4),b(2,4),b(3,4))
      return
14    engclc=elj(rr*(cutof2/cut2))
      return
15    engclc=eijljg(c(1,i1),b)
      return
16    engclc=eijgen(c(1,i1),b)
      return
      end
      function ewwmcy(
     -  a1x,a1y,a1z,a2x,a2y,a2z,a3x,a3y,a3z,a4x,a4y,a4z,
     -  b1x,b1y,b1z,b2x,b2y,b2z,b3x,b3y,b3z,b4x,b4y,b4z)
c*****Calculates ab initio ci interaction energy between two
c*****water molecules.  see O.Matsuoka, E.Clementi and
c*****M.Yoshimine,  JCP 64,1351-1361(1976).
c     The atomic locations in the article correspond to
c     h1=i+1, h2=i+2, h3=j+1, h4=j+2, o5=i, o6=j.
      r562x=a1x-b1x
      r562y=a1y-b1y
      r562z=a1z-b1z
      r562=r562x*r562x+r562y*r562y+r562z*r562z
      r56=sqrt(r562)
      r562x=a2x-b2x
      r562y=a2y-b2y
      r562z=a2z-b2z
      r562=r562x*r562x+r562y*r562y+r562z*r562z
      r13=sqrt(r562)
      r562x=a2x-b3x
      r562y=a2y-b3y
      r562z=a2z-b3z
      r562=r562x*r562x+r562y*r562y+r562z*r562z
      r14=sqrt(r562)
      r562x=a3x-b2x
      r562y=a3y-b2y
      r562z=a3z-b2z
      r562=r562x*r562x+r562y*r562y+r562z*r562z
      r23=sqrt(r562)
      r562x=a3x-b3x
      r562y=a3y-b3y
      r562z=a3z-b3z
      r562=r562x*r562x+r562y*r562y+r562z*r562z
      r24=sqrt(r562)
      r562x=a2x-b4x
      r562y=a2y-b4y
      r562z=a2z-b4z
      r562=r562x*r562x+r562y*r562y+r562z*r562z
      r18=sqrt(r562)
      r562x=a3x-b4x
      r562y=a3y-b4y
      r562z=a3z-b4z
      r562=r562x*r562x+r562y*r562y+r562z*r562z
      r28=sqrt(r562)
      r562x=a4x-b2x
      r562y=a4y-b2y
      r562z=a4z-b2z
      r562=r562x*r562x+r562y*r562y+r562z*r562z
      r37=sqrt(r562)
      r562x=a4x-b3x
      r562y=a4y-b3y
      r562z=a4z-b3z
      r562=r562x*r562x+r562y*r562y+r562z*r562z
      r47=sqrt(r562)
      r562x=a4x-b4x
      r562y=a4y-b4y
      r562z=a4z-b4z
      r562=r562x*r562x+r562y*r562y+r562z*r562z
      r78=sqrt(r562)
c     Calculate interaction energy.
      ewwmcy=+qq*(1.0/r13+1.0/r14+1.0/r23+1.0/r24+4.0/r78)
      ewwmcy=ewwmcy -dqq*((1.0/r18+1.0/r28+1.0/r37+1.0/r47))
c     Calculate remaining distances.
      if (r56 .ge. 11.0) return
      r562x=a2x-b1x
      r562y=a2y-b1y
      r562z=a2z-b1z
      r562=r562x*r562x+r562y*r562y+r562z*r562z
      r16=sqrt(r562)
      r562x=a3x-b1x
      r562y=a3y-b1y
      r562z=a3z-b1z
      r562=r562x*r562x+r562y*r562y+r562z*r562z
      r26=sqrt(r562)
      r562x=a1x-b2x
      r562y=a1y-b2y
      r562z=a1z-b2z
      r562=r562x*r562x+r562y*r562y+r562z*r562z
      r35=sqrt(r562)
      r562x=a1x-b3x
      r562y=a1y-b3y
      r562z=a1z-b3z
      r562=r562x*r562x+r562y*r562y+r562z*r562z
      r45=sqrt(r562)
      ewwmcy=ewwmcy+(a1*exp(-b1*r56))
      ewwmcy=ewwmcy +(a2*(exp(-b2*r13)+exp(-b2*r14)
     -         +exp(-b2*r23)+exp(-b2*r24)))
      ewwmcy=ewwmcy +(a3*(exp(-b3*r16)+exp(-b3*r26)
     -         +exp(-b3*r35)+exp(-b3*r45)))
      ewwmcy=ewwmcy-(a4*(exp(-b4*r16)+exp(-b4*r26)
     -        +exp(-b4*r35)+exp(-b4*r45)))
      return
      end
      function ewwst2(a,b)
c*****ST2 water-water cyclic coding
      dimension a(3,5),b(3,5)
      common /st2par/ qww(5,5),rl,rl2,ru2,r3ul,rul3i,sg6,eps4s
c     ru2=ru**2; r3ul=3*ru-rl; sg6=sigma**6; eps4s=4*eps*sigma**6
      common /info/ nslv,nmolg,nmolr,natg,natr,icopvg,icopvr,nmolmx
      real*8 vq
c     Atoms are in the order : o, h1, h2, lp1, lp2.
      roo2=(a(1,1)-b(1,1))**2+(a(2,1)-b(2,1))**2+(a(3,1)-b(3,1))**2
      r6=1.0/(roo2*roo2*roo2)
      r12=r6*r6
      ewwst2=eps4s*(sg6*r12-r6)
      if (roo2 .lt. rl2) return
      s=1.0
      if (roo2 .gt. ru2) go to 15
      roo=sqrt(roo2)
      s=(roo-rl)**2*(r3ul-roo-roo)*rul3i
15    vq=0.d0
      do 20 i=2,nslv
        do 20 j=2,nslv
20        vq=vq+dble(qww(i,j)/sqrt((a(1,i)-b(1,j))**2+
     -      (a(2,i)-b(2,j))**2+(a(3,i)-b(3,j))**2))
      vqs=vq
      ewwst2=ewwst2+s*vqs
      return
      end
      function eijljg(a,b)
c*****General 6-12 molecule, cyclic coding
      dimension a(3,10),b(3,10)
      common /info/ nslv,nmolg,nmolr,natg,natr,icopvg,icopvr,nmolmx
      common /genslv/ c6(10),c12(10),cc6(10,10),cc12(10,10),qqg(10,10),
     -                siglj(10),epslj(10)
      real*8 vq
15    vq=0.d0
      do 20 i=1,nslv
        do 20 j=1,nslv
          rr=(a(1,i)-b(1,j))**2+
     -      (a(2,i)-b(2,j))**2+(a(3,i)-b(3,j))**2
          r6=rr*rr*rr
          r12=r6*r6
20        vq=vq+cc12(i,j)/r12-cc6(i,j)/r6
      eijljg=vq
      return
      end
      function eijgen(a,b)
c*****General solvent cyclic coding
      dimension a(3,10),b(3,10)
      character*80 label
      common /id/ label(2),temp,fkt,cutof,cutof2,edgea,edge,hlfedg,
     -  prsag,prsar,pressg,pressr,bparg,bparr
      common /info/ nslv,nmolg,nmolr,natg,natr,icopvg,icopvr,nmolmx
      common /genslv/ c6(10),c12(10),cc6(10,10),cc12(10,10),qqg(10,10),
     -                siglj(10),epslj(10)
      common /option/ iop(40)
      real*8 vq
15    vq=0.d0
      do 20 i=1,nslv
        do 20 j=1,nslv
          rr=(a(1,i)-b(1,j))**2+
     -      (a(2,i)-b(2,j))**2+(a(3,i)-b(3,j))**2
          r6=rr*rr*rr
          r12=r6*r6
          if(iop(9).eq.0) then
            vq=vq+dble(qqg(i,j)/sqrt(rr))+cc12(i,j)/r12-cc6(i,j)/r6
          else
            rr1=sqrt(rr)
            cfac=1.+0.5*(rr1/cutof)**3
            vq=vq+dble(cfac*qqg(i,j)/rr1)+cc12(i,j)/r12-cc6(i,j)/r6      
          endif
20    continue
      eijgen=vq
      return
      end
      function ewwspc(
     -  a1x,a1y,a1z,a2x,a2y,a2z,a3x,a3y,a3z,
     -  b1x,b1y,b1z,b2x,b2y,b2z,b3x,b3y,b3z)
c*****SPC (or TIP3P) water
      common /spcpar/qqspc,qqspc2,qqspc4,alj,blj,alj6,alj66,blj12,blj122
      r562x=a1x-b1x
      r562y=a1y-b1y
      r562z=a1z-b1z
      r562=r562x*r562x+r562y*r562y+r562z*r562z
      r56=sqrt(r562)
      r566=r562*r562*r562
      r5612=r566*r566
      r562x=a2x-b2x
      r562y=a2y-b2y
      r562z=a2z-b2z
      r562=r562x*r562x+r562y*r562y+r562z*r562z
      r13=sqrt(r562)
      r562x=a2x-b3x
      r562y=a2y-b3y
      r562z=a2z-b3z
      r562=r562x*r562x+r562y*r562y+r562z*r562z
      r14=sqrt(r562)
      r562x=a3x-b2x
      r562y=a3y-b2y
      r562z=a3z-b2z
      r562=r562x*r562x+r562y*r562y+r562z*r562z
      r23=sqrt(r562)
      r562x=a3x-b3x
      r562y=a3y-b3y
      r562z=a3z-b3z
      r562=r562x*r562x+r562y*r562y+r562z*r562z
      r24=sqrt(r562)
      r562x=a2x-b1x
      r562y=a2y-b1y
      r562z=a2z-b1z
      r562=r562x*r562x+r562y*r562y+r562z*r562z
      r16=sqrt(r562)
      r562x=a3x-b1x
      r562y=a3y-b1y
      r562z=a3z-b1z
      r562=r562x*r562x+r562y*r562y+r562z*r562z
      r26=sqrt(r562)
      r562x=a1x-b2x
      r562y=a1y-b2y
      r562z=a1z-b2z
      r562=r562x*r562x+r562y*r562y+r562z*r562z
      r53=sqrt(r562)
      r562x=a1x-b3x
      r562y=a1y-b3y
      r562z=a1z-b3z
      r562=r562x*r562x+r562y*r562y+r562z*r562z
      r54=sqrt(r562)
c     Calculate interaction energy.
      ewwspc=(qqspc/r13+qqspc/r14+qqspc/r23+qqspc/r24+qqspc4/r56)
      ewwspc=ewwspc -(qqspc2/r16+qqspc2/r26+qqspc2/r53+qqspc2/r54)
      ewwspc=ewwspc-alj6/r566
      ewwspc=ewwspc+blj12/r5612
      return
      end
      function ewwtip(
     -  a4x,a4y,a4z,a2x,a2y,a2z,a3x,a3y,a3z,a1x,a1y,a1z,
     -  b4x,b4y,b4z,b2x,b2y,b2z,b3x,b3y,b3z,b1x,b1y,b1z)
c*****TIP4P(2p) water
      common /spcpar/qqtip,qqtip2,qqtip4,alj,blj,alj6,alj66,blj12,blj122
      r562x=a1x-b1x
      r562y=a1y-b1y
      r562z=a1z-b1z
      r562=r562x*r562x+r562y*r562y+r562z*r562z
      r56=sqrt(r562)
      r562x=a4x-b4x
      r562y=a4y-b4y
      r562z=a4z-b4z
      r562=r562x*r562x+r562y*r562y+r562z*r562z
      r566=r562*r562*r562
      r5612=r566*r566
      r562x=a2x-b2x
      r562y=a2y-b2y
      r562z=a2z-b2z
      r562=r562x*r562x+r562y*r562y+r562z*r562z
      r13=sqrt(r562)
      r562x=a2x-b3x
      r562y=a2y-b3y
      r562z=a2z-b3z
      r562=r562x*r562x+r562y*r562y+r562z*r562z
      r14=sqrt(r562)
      r562x=a3x-b2x
      r562y=a3y-b2y
      r562z=a3z-b2z
      r562=r562x*r562x+r562y*r562y+r562z*r562z
      r23=sqrt(r562)
      r562x=a3x-b3x
      r562y=a3y-b3y
      r562z=a3z-b3z
      r562=r562x*r562x+r562y*r562y+r562z*r562z
      r24=sqrt(r562)
      r562x=a2x-b1x
      r562y=a2y-b1y
      r562z=a2z-b1z
      r562=r562x*r562x+r562y*r562y+r562z*r562z
      r16=sqrt(r562)
      r562x=a3x-b1x
      r562y=a3y-b1y
      r562z=a3z-b1z
      r562=r562x*r562x+r562y*r562y+r562z*r562z
      r26=sqrt(r562)
      r562x=a1x-b2x
      r562y=a1y-b2y
      r562z=a1z-b2z
      r562=r562x*r562x+r562y*r562y+r562z*r562z
      r53=sqrt(r562)
      r562x=a1x-b3x
      r562y=a1y-b3y
      r562z=a1z-b3z
      r562=r562x*r562x+r562y*r562y+r562z*r562z
      r54=sqrt(r562)
c     Calculate interaction energy.
      ewwtip=(qqtip/r13+qqtip/r14+qqtip/r23+qqtip/r24+qqtip4/r56)
      ewwtip=ewwtip -(qqtip2/r16+qqtip2/r26+qqtip2/r53+qqtip2/r54)
      ewwtip=ewwtip-alj6/r566
      ewwtip=ewwtip+blj12/r5612
      return
      end
      function elj(rr)
c*****lennard-rJnes
      rr6=rr*rr*rr
      r6=4.0/rr6
      elj=r6/rr6-r6
      return
      end
      subroutine mnmxup(x,nmc,xmin,xmax,nxmin,nxmax)
c*****Updates the extrema values and the corresponding stepnumbers
      if (x .lt. xmin) then
        xmin=x
        nxmin=nmc
      else if (x .gt. xmax) then
        xmax=x
        nxmax=nmc
      end if
      return
      end
      subroutine trnsfi(l,m,n)
c*****Fast transfer - integer arrays
      dimension l(n),m(n)
      do i=1,n
        l(i)=m(i)
      end do
      return
      end
      subroutine trnsfr(a,b,n)
c*****Fast transfer
      dimension a(n),b(n)
      do i=1,n
        a(i)=b(i)
      end do
      return
      end
      subroutine vprd(a,b,c)
c*****Computes the vector product a x b and saves it into c
      dimension a(3),b(3),c(3)
      c(1)=a(2)*b(3)-b(2)*a(3)
      c(2)=a(3)*b(1)-b(3)*a(1)
      c(3)=a(1)*b(2)-b(1)*a(2)
      return
      end
      subroutine vprod(r,i,j,k)
c*****Computes the vector product of the columns i and j into k
      dimension r(3,3)
      do l=1,3
        r(l,k)=r(mod(l,3)+1,i)*r(mod(l+1,3)+1,j)-
     -         r(mod(l,3)+1,j)*r(mod(l+1,3)+1,i)
      end do
      return
      end
      subroutine mnorm(r)
c*****Normalizes the matrix r
      dimension r(3,3)
      do i=1,3
        rr=0.0
        do k=1,3
          rr=rr+r(k,i)**2
        end do
        rr=sqrt(rr)
        do k=1,3
          r(k,i)=r(k,i)/rr
        end do
      end do
      return
      end
      subroutine disfin
c*****Compute the distribution functions, save the plot file
      common /coord/ relcrd(3,10),cslv(3,10),qslv(10),islv(10)
      real*8 ug,ur,uug,uur,vg,vr,vvg,vvr,pg,pr,ppg,ppr,dva,dvag,dvar
      common /accu/ ug,ur,uug,uur,vg,vr,vvg,vvr,pg,pr,ppg,ppr,dva,
     -  dvag,dvar
      common /svint/ nmc,nmaccg,nmaccr,neming,neminr,nemaxg,nemaxr,
     -  nmcmax,nmcrep,nrecd,nmcvch,nsync,niaccg,ndaccg,nmovg,nmovr,
     -  nvming,nvmaxg,npming,npmaxg,nvminr,nvmaxr,npminr,npmaxr,
     -  nvacc,nvaccg,nvaccr,nmccnt,ncnt
      common /sreal/ eming,eminr,emaxg,emaxr,vol,volg,volr,vfacg,vfacr,
     -  cut2g,cut2r,vming,vmaxg,vminr,vmaxr
      common /dploy/ rsi(3,10),crmi(3),orj(3,3),cedgs,rotax,delvol
      common /info/ nslv,nmolg,nmolr,natg,natr,icopvg,icopvr,nmolmx
      common /nmbrs/ autokc,autoa,rdtodg,boltz,pi,pi2
      common /atwght/ avgd,aw(100)
      character*80 label
      common /id/ label(2),temp,fkt,cutof,cutof2,edgea,edge,hlfedg,
     -  prsag,prsar,pressg,pressr,bparg,bparr
      common /grid/ ncover(100,100,100),igfree(1000000),ngfree,ngrid,
     -  ng1,ng2,ng3,fltngr,fltng3,grid,e2g,e2gm,rsppe2,rspme2
      common /cavdat/ rsph,rsphs,corlng,edgrs,edgrs2,edgrq2
      common /cavprb/ ndistg(1000),ndistr(1000),pming,pmaxg,pminr,pmaxr
      real*8 angcrg,angcrr
      common /distr/ angcrg(210),angcrr(210),nrg(210),nrr(210),
     -  rdbg(10,210),rdbr(10,210),angcbg(10,210),angcbr(10,210),
     -  ncnfbg(10),ncnfbr(10),rgrid,rgridi,rgrdig,rgrdir,rdmax2,rdmx2g,
     -  rdmx2r,rdm2ng,rdm2nr,rgring,rgrinr,nrgrid,icolcr,nconfg,nconfr,
     -  lenblk,nblk
      common /option/ iop(40)
      common /cntrl/ ceg(100),cer(100),cvg(100),cng(100)
      character*8 itg,itr
      common /label/ itg,itr
      character*8 ibl
      character*12 iang,irdf
      dimension uc(100)
      data ibl/'        '/
      data iang/'ang. correl. '/,irdf/'radial dist. '/
      write (6,1000) nmc,id,temp
      if (nmc .eq. 0) return
      call csave(10)
      nmcv=nmc/nmcvch
      nmce=nmc+nmcv
      if (iop(6) .gt. 0) then
        nmce=nmce+nmc
        pavg=pg/nmc
        pavr=pr/nmc
      else
        pavg=nmolg
        pavr=nmolr
      end if
      uavg=ug/nmce
      usdg=dsqrt(dabs(uug/nmce-uavg**2))
      uavr=ur/nmce
      usdr=dsqrt(dabs(uur/nmce-uavr**2))
      ug=uavg/pavg
      ur=uavr/pavr
      write (6,1001) itg,uavg,usdg,eming,neming,emaxg,nemaxg,ug
      write (6,1001) itr,uavr,usdr,eminr,neminr,emaxr,nemaxr,ur
c     Calculate error estimates
      nbl=nmccnt+nmccnt/nmcvch
      if (iop(6) .gt. 0) nbl=nbl+nmccnt
      call blckav(ceg,uc,nbl,ncnt)
      call btchmn(ncnt,nmccnt,uc,' g: <e>=',6,0,sd2)
      if (pmaxr .gt. 0) then
        call blckav(cer,uc,nbl,ncnt)
        call btchmn(ncnt,nmccnt,uc,' r: <e>=',6,0,sd2)
      end if
      if (nmc .ge. nmcvch) then
        vavg=vg/nmcv
        vsdg=dsqrt(dabs(vvg/nmcv-vavg**2))
        vavr=vr/nmcv
        vsdr=dsqrt(dabs(vvr/nmcv-vavr**2))
        write (6,1003) itg,vavg,vsdg,vming,nvming,vmaxg,nvmaxg
        write (6,1003) itr,vavr,vsdr,vminr,nvminr,vmaxr,nvmaxr
c       Calculate error estimates
        nbl=nmccnt/nmcvch
        call blckav(cvg,uc,nbl,ncnt)
        call btchmn(ncnt,nmccnt,uc,' g: <v>=',6,0,sd2v)
        if (iop(3) .eq. 0) then
          pacc=float(nvacc)/float(nmcv)
          dvavg=0
          if (nvacc .gt. 0) dvavg=dva/nvacc
          write (6,1004) ibl,dvavg,pacc
        else if (iop(3) .eq. 1) then
          paccg=float(nvaccg)/float(nmcv)
          dvavg=0
          if (nvaccg .gt. 0) dvavg=dvag/nvaccg
          paccr=float(nvaccr)/float(nmcv)
          dvavr=0
          if (nvaccr .gt. 0) dvavr=dvar/nvaccr
          write (6,1004) itg,dvavg,paccg,itr,dvavr,paccr
        end if
      end if
      if (iop(6) .gt. 0) then
        psdg=dsqrt(dabs(ppg/nmc-pavg**2))
        psdr=dsqrt(dabs(ppr/nmc-pavr**2))
        nming=pming
        nmaxg=pmaxg
        nminr=pminr
        nmaxr=pmaxr
        write (6,1005) itg,pavg,psdg,nming,npming,nmaxg,npmaxg
        write (6,1005) itr,pavr,psdr,nminr,npminr,nmaxr,npmaxr
c       Calculate error estimates
        nbl=nmccnt
        call blckav(cng,uc,nbl,ncnt)
        call btchmn(ncnt,nmccnt,uc,' g: <n>=',6,0,sd2n)
        write (6,1007) itg,nming,nmaxg,(ndistg(i),i=nming,nmaxg)
        if (nminr .gt. 0) then
          write (6,1007) itr,nminr,nmaxr,(ndistr(i),i=nminr,nmaxr)
        else
c         nmolr=0 was not accumulated
          write (6,1007) itr,nminr,nmaxr,
     -      (ndistg(nmaxg-i+nming),i=nming,nmaxg)
        end if
        pacc=float(niaccg+ndaccg)/float(nmc)
        write (6,1006) pacc,niaccg,ndaccg
        if (iop(6) .gt. 1) then
          pcav=ngfree/fltng3
          write (6,1008) rsph,pcav,ngrid,ngfree
        end if
      end if
      if (nmcvch .gt. nmc) then
        vavg=volg
        vavr=volr
      end if
      if (iop(8) .eq. 1) then
        corlg=0.0
        corlr=0.0
        if (vavg .gt. 0.0) corlg=corlng*pavg**2/vavg
        if (vavr .gt. 0.0) corlr=corlng*pavr**2/vavr
        corl1g=corlg/pavg
        ulg=(uavg+corlg)/pavg
        corl1r=corlr/pavr
        ulr=(uavr+corlr)/pavr
        write (6,1009) itg,corlg,corl1g,ulg
        write (6,1009) itr,corlr,corl1r,ulr
      end if
      deng=pavg/vavg
      denr=pavr/vavr
      releg=0.0
      if (pavg .gt. 0.0) releg=releg+sd2n/pavg
      if (vavg .gt. 0.0) releg=releg+sd2v/vavg
      reler=0.0
      if (pavr .gt. 0.0) reler=reler+sd2n/pavr
      if (vavr .gt. 0.0) reler=reler+sd2v/vavr
      denge=deng*releg
      denre=denr*reler
      write (6,1010) itg,itr,deng,denge,denr,denre
      if (iop(11) .ne. 4) then
        rmw=0.0
        do 30 k=1,nslv
30        rmw=rmw+aw(islv(k))
        denfac=rmw/avgd/1.e-24
        dengg=deng*denfac
        denrr=denr*denfac
        denge=denge*denfac
        denre=denre*denfac
        write (6,1011) itg,itr,dengg,denge,denrr,denre,rmw
      end if
c     Estimating n or v for equivolume solution
      ng1=vol*deng
      nr1=vol*denr
      voln=vol*float(nmolg+nmolr)/float(ng1+nr1)
      edgean=voln**(1.0/3.0)
      write (6,1012) ng1,nr1,voln,edgean
      if (iop(17) .gt. 0) then
c       Print rdf results 
        write (6,1015) nconfg,itg,nconfr,itr
        if (iop(17) .gt. 1) then
          do i=1,nrgrid
            if (nrg(i) .gt. 0) angcrg(i)=angcrg(i)/nrg(i) 
            if (nblk .gt. 0) then
              do ib=1,nblk
                if (rdbg(ib,i) .ne. 0.0)
     -            angcbg(ib,i)=ncnfbg(ib)*angcbg(ib,i)/rdbg(ib,i)
              end do
            end if
          end do
          write (6,1013) itg,(angcrg(i),i=1,nrgrid)
          if (nconfr .gt. 0) then
            do i=1,nrgrid
              if (nrr(i) .gt. 0) angcrr(i)=angcrr(i)/nrr(i)  
              if (nblk .gt. 0) then
                do ib=1,nblk
                  if (rdbr(ib,i) .ne. 0.0)
     -              angcbr(ib,i)=ncnfbr(ib)*angcbr(ib,i)/rdbr(ib,i)
                end do
              end if
            end do
            write (6,1013) itr,(angcrr(i),i=1,nrgrid)
          end if
        end if
        do i=1,nrgrid
          rr=(i-1)*rgrid+rgrid/2.0
          gridfac=2.0/(deng*4.0*pi*rr**2*rgrid*pavg)
          angcrg(i)=nrg(i)*gridfac/nconfg
          if (nblk .gt. 0) then
            do ib=1,nblk
             rdbg(ib,i)=rdbg(ib,i)*gridfac
            end do
          end if
        end do
        write (6,1014) itg,(angcrg(i),i=1,nrgrid)
        if (nconfr .gt. 0) then
          do i=1,nrgrid
            rr=(i-1)*rgrid+rgrid/2.0
            gridfac=2.0/(denr*4.0*pi*rr**2*rgrid*pavr)
            angcrr(i)=nrr(i)*gridfac/nconfr
            if (nblk .gt. 0) then
              do ib=1,nblk
               rdbr(ib,i)=rdbr(ib,i)*gridfac
              end do
            end if
          end do
          write (6,1014) itr,(angcrr(i),i=1,nrgrid)
        end if
      end if
      if (nblk .gt. 1) then
c       Calculate error estimates
        call rdferr(nrgrid,lenblk,nblk,itg,irdf,rdbg,ncnfbg)
        call rdferr(nrgrid,lenblk,nblk,itr,irdf,rdbr,ncnfbr)
        if (iop(17) .gt. 1) then
          call rdferr(nrgrid,lenblk,nblk,itg,iang,angcbg,ncnfbg)
          call rdferr(nrgrid,lenblk,nblk,itr,iang,angcbr,ncnfbr)
        end if
      end if
      call crstr(10)
      return
1000  format(/,' Results at MC step number ',i8,':',2(/,1x,20a4),/,
     -  ' temperature=',f8.2,' Kelvin')
1001  format(1x,a8,' <e>=',e12.5,' sd=',e12.5,' min=',e12.5,' (',
     -  i8,') max=',e12.5,'(',i8,') <e>/<n>=',f12.5)
1003  format(1x,a8,' <v>=',f12.4,' sd=',f12.4,' min=',f12.4,' (',
     -  i8,') max=',f12.4,'(',i8,')')
1004  format(1x,a8,' average volume change size=',f8.4,' A**3',
     -  '   Volume exchange acceptance rate=',f6.3)
1005  format(1x,a8,' <n>=',f12.5,' sd=',f12.5,' min=',i12,' (',
     -  i8,') max=',i12,'(',i8,')')
1006  format(' Particle exchange acceptance rate=',f12.8,
     -  '  number of accepted insertions and deletions=',2i9)
1007  format(1x,a8,' the distribution of the number of molecules ',
     -  'from ',i4,' to ',i4,':',/,(15i8))
1008  format(' Cavity-biased insertions, cavity radius=',f6.2,' A ',
     -  ' probability of finding a cavity=',f8.6,' ngrid=',i3,
     -  ' ngfree=',i6)
1009  format(1x,a8,' long-range corrections (total)=',e12.5,
     -  ' (per molecule)=',f10.5,' corrected energy=',e12.5)
1010  format(/,' The average densities in the ',a8,' and ',a8,' :',
     -  2(f10.6,' +/-',f8.6),' n/A**3')
1011  format(/,' The average densities in the ',a8,' and ',a8, ':',
     -  2(f10.6,' +/-',f8.6),' gram/ml  molecular weight=',f9.4)
1012  format(' To obtain equivolume systems, use ',i4,' and ',i4,
     -  ' particles into systems g and r, respectively, or,',/,
     -  ' change the initial common volume to ',f10.4,' A**3.',
     -  ' this would require edge=',f12.6,' A')
1013  format(/,' The average angular correlations as a function of',
     -  ' the distance in the ',2a4,':',/,(20f6.3))
1014  format(/,' The radial distribution function of the ',a8,
     -  ':',/,(20f6.2))
1015  format(/,' Distribution functions were based on ',i6,
     -  ' configurations in the ',a8,' and ',i6,
     -  ' configurations in the ',a8)
      end
      subroutine rdferr(nrgrid,lenblk,nblk,isyst,itype,cumsum,ncnfb)
      dimension cumsum(10,210),sd(210),itype(3),ncnfb(10),blav(10)
      character*12 itype
      character*8 isyst
c*****Calculate error estimate for distribution functions
      do i=1,nrgrid
        call blckav(cumsum(1,i),blav,lenblk,nblk)
c       Correct for possible uneven blocklength
        do ib=1,nblk
          len=ncnfb(ib)
          if (ib .gt. 1) len=len-ncnfb(ib-1)
          if (len .gt. 0) blav(ib)=blav(ib)*float(lenblk)/float(len)
        end do
        call btchmn(nblk,lenblk,blav,isyst,6,1,sd(i))
      end do
      write (6,1000) itype,isyst,(sd(i),i=1,nrgrid)
      return
1000  format(/,' The error estimates on ',a12,' for the ',a8,' :',/,
     -  (20f6.3))
      end
      subroutine crorgn(edge)
c*****Initialize the image cells
      common /ecell/cic(3,27)
      dimension y(3),x(3)
c     cic contains the negative of the respective cell center coords
c     siple cubic
      y(1)=0.
      y(2)=edge
      y(3)=- edge
      ii=0
      do i=1,3
        x(3)=y(i)
        do k=1,3
          x(2)=y(k)
          do l=1,3
            x(1)=y(l)
            ii=ii+1
            do j=1,3
              cic(j,ii)=x(j)
            end do
          end do
        end do
      end do
      return
      end
      subroutine cellfn(y1,y2,y3,x1,x2,x3,hlfedg,ipbc)
c*****Find the cell where the point (x1,x2,x3) is and generate the
c*****displacement moving it back to the first cell.
      common /option/ iop(40)
      common /ecell/cic(81)
c     Simple cubic
      hlfedm=-hlfedg
      n1=0
      if (x1 .ge. hlfedg) n1=2
      if (x1. le. hlfedm) n1=1
      n2=0
      if (x2 .ge. hlfedg) n2=2
      if (x2 .le. hlfedm) n2=1
      n3=0
      if (x3 .ge. hlfedg) n3=2
      if (x3 .le. hlfedm) n3=1
      ipbc=n1+n2*3+n3*9
      jdim=ipbc+ipbc+ipbc
      ipbc=ipbc+1
      y1=cic(jdim+1)
      y2=cic(jdim+2)
      y3=cic(jdim+3)
      return
      end
      subroutine movex(texo,ex,isol,nmol)
c*****Move the nmol-th molecules pair energies to the place of the imolth
      dimension texo(1000),ex(499500)
      nmol1=nmol-1
      isol0=(isol-2)*(isol-1)/2
      do k=1,nmol1
        texo(k)=ex(isol0+k)
      end do
      do k=1,isol-1
        ex(isol0+k)=texo(k)
      end do
      isol1=(isol-1)*(isol)/2+isol
      do k=isol+1,nmol
        ex(isol1)=texo(k)
        isol1=isol1+k-1
      end do
      return
      end
      subroutine savegy(tesi,tex,texo,ex,isol,nmol)
c*****Put the pair energies into ex and update the binding energy
c*****array, tesi
      real*8 tesi,corect
      dimension tesi(1000),tex(1000),texo(1000),ex(499500)
      isol0=(isol-2)*(isol-1)/2
      do k=1,isol-1
        corect= dble(tex(k))-dble(texo(k))
        tesi(k)=tesi(k)+corect
        tesi(isol)=tesi(isol)+corect
        ex(isol0+k)=tex(k)
      end do
      isol1=(isol-1)*(isol)/2+isol
      do k=isol+1,nmol
        corect= dble(tex(k))-dble(texo(k))
        tesi(k)=tesi(k)+corect
        tesi(isol)=tesi(isol)+corect
        ex(isol1)=tex(k)
        isol1=isol1+k-1
      end do
      return
      end
      subroutine getegy(texo,ex,isol,nmol)
c*****Put eww(isol,i) into texo(i) from the saved pair energies (ex)
      dimension texo(1000),ex(499500)
c     If the number of molecules is greater than, say, 300 then
c     ex array is to be eliminated and this routine will have to
c     call eww instead.
      texo(isol)=0.0
      isol0=(isol-2)*(isol-1)/2
      do k=1,isol-1
        texo(k)=ex(isol0+k)
      end do
      isol1=(isol-1)*(isol)/2+isol
      do k=isol+1,nmol
        texo(k)=ex(isol1)
        isol1=isol1+k-1
      end do
      return
      end
      subroutine blckav(cumsum,blav,nblk,n)
c*****Converts cumulative sums to block averages
      dimension cumsum(n),blav(n)
      if (n .lt. 1) return
      blav(1)=cumsum(1)
      if (n .gt. 1) then
        do i=2,n
          blav(i)=cumsum(i)-cumsum(i-1)
      end do
      end if
      do i=1,n
        blav(i)=blav(i)/nblk
      end do
      return
      end
      subroutine btchmn(npts,nintvl,uinp,inam,iout,iprt,sd2i)
c#    MMC routine 082 lstmod: 05/09/91
c*****Computes error bound with the method of batch means
      dimension uinp(npts),u(100)
      character*8 inam
      dimension nmncrt(20,20),nmxcrt(20,20)
      character*12 uncorr,corr,high,low,decide
      data uncorr/'uncorrelated'/,corr/'correlated  '/,
     -  high/' >>>        '/,low/' ???        '/
c     Minimum critical values
      data nmncrt/20*0,
     -  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2,
     -  0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3,
     -  0, 0, 0, 0, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4,
     -  0, 0, 0, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5,
     -  0, 0, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6,
     -  0, 0, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6,
     -  0, 0, 2, 3, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 7, 7,
     -  0, 0, 2, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8,
     -  0, 0, 2, 3, 3, 4, 5, 5, 5, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 9,
     -  0, 0, 2, 3, 4, 4, 5, 5, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 9,
     -  0, 2, 2, 3, 4, 4, 5, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9,10,10,
     -  0, 2, 2, 3, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 9,10,10,10,10,
     -  0, 2, 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 9, 9,10,10,10,11,11,
     -  0, 2, 3, 3, 4, 5, 6, 6, 7, 7, 8, 8, 9, 9,10,10,11,11,11,12,
     -  0, 2, 3, 4, 4, 5, 6, 6, 7, 8, 8, 9, 9,10,10,11,11,11,12,12,
     -  0, 2, 3, 4, 4, 5, 6, 7, 7, 8, 9, 9,10,10,11,11,11,12,12,13,
     -  0, 2, 3, 4, 5, 5, 6, 7, 8, 8, 9, 9,10,10,11,11,12,12,13,13,
     -  0, 2, 3, 4, 5, 6, 6, 7, 8, 8, 9,10,10,11,11,12,12,13,13,13,
     -  0, 2, 3, 4, 5, 6, 6, 7, 8, 9, 9,10,10,11,12,12,13,13,13,14/
c     Maximum critical values:
      data nmxcrt/60*0,
     -  0, 0, 0, 0, 9, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     -  0, 0, 0, 9,10,10,11,11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     -  0, 0, 0, 9,10,11,12,12,13,13,13,13, 0, 0, 0, 0, 0, 0, 0, 0,
     -  0, 0, 0, 0,11,12,12,12,14,14,14,14,15,15,15, 0, 0, 0, 0, 0,
     -  0, 0, 0, 0,11,12,13,13,14,15,15,16,16,16,16,17,17,17,17,17,
     -  0, 0, 0, 0, 0,13,14,14,15,16,16,16,17,17,18,18,18,18,18,18,
     -  0, 0, 0, 0, 0,13,14,15,16,16,17,17,18,18,18,19,19,19,20,20,
     -  0, 0, 0, 0, 0,13,14,15,16,17,17,18,19,19,19,20,20,20,21,21,
     -  0, 0, 0, 0, 0,13,14,16,16,17,18,19,19,20,20,21,21,21,22,22,
     -  0, 0, 0, 0, 0, 0,15,16,17,18,19,19,20,20,21,21,22,22,23,23,
     -  0, 0, 0, 0, 0, 0,15,16,17,18,19,20,20,21,22,23,23,23,23,24,
     -  0, 0, 0, 0, 0, 0,15,16,18,18,19,20,21,22,22,23,23,24,24,25,
     -  0, 0, 0, 0, 0, 0, 0,17,18,19,29,21,21,22,23,23,24,25,25,25,
     -  0, 0, 0, 0, 0, 0, 0,17,18,19,20,21,22,23,23,24,25,25,26,26,
     -  0, 0, 0, 0, 0, 0, 0,17,18,19,20,21,22,23,24,25,25,26,26,27,
     -  0, 0, 0, 0, 0, 0, 0,17,18,29,21,22,23,23,24,25,26,26,27,27,
     -  0, 0, 0, 0, 0, 0, 0,17,18,20,21,22,23,25,25,26,26,27,27,28/
      sd2i=0.0
      if (npts .lt. 2) return
      call trnsfr(u,uinp,npts)
      npts0=npts
      nbl=nintvl
      do il=1,10
        call var(u,npts0,sd2,nup,ndown,nrun)
        if (il .eq. 1) sd2i=sd2
        if (ndown .gt. 20 .or. nup .gt. 20) then
          decide=high
        else if (ndown .eq. 0 .or. nup .eq. 0) then
          decide=low
        else if (nmxcrt(ndown,nup) .eq. 0 .or.
     -           nmncrt(ndown,nup) .eq. 0) then
          decide=low
        else if (nrun .lt. nmncrt(ndown,nup) .or.
     -           nrun .gt. nmxcrt(ndown,nup)) then
          decide=corr
        else
          decide=uncorr
        end if
        if (iprt .eq. 0) write (iout,1002) inam,sd2,nbl,nup,ndown,
     -                   nrun,decide
c       Double the blocksize
        nbl=nbl*2
        npts0=npts0/2
        if (npts0 .le. 2) return
c       Clump together every pair of subsequent block averages
        do j=1,npts0
          u(j)=(u(2*j)+u(2*j-1))/2.0
        end do
      end do
      return
1002  format(1x,a8,' 2*SD=',f10.4,' block size=',i8,
     -  ' nup=',i3,' ndown=',i3,' nrun=',i3,2x,a12)
      end
      subroutine var(a,n,sd2,nup,ndown,nrun)
c#    MMC routine 083 lstmod: 10/22/85
c*****Computes the variance for the method of batch means
      dimension a(n),rs(100)
      real*8 av,sd
      sd=0.0d0
      av=0.0d0
      do 2 i=1,n
        av=av+a(i)
2       sd=sd+a(i)**2
      sd=sqrt(abs(sngl(sd)/n-(sngl(av)/n)**2)/float(n-1))
      sd2=2.0*sd
c     Find the median
      do 23 i=1,n
23      rs(i)=a(i)
      n2=n/2+1
      do 20 i=1,n2
        r=rs(i)
        jm=i
        do 21 j=i,n
          if (rs(j) .le. r) go to 21
          r=rs(j)
          jm=j
21      continue
        rs(jm)=rs(i)
20      rs(i)=r
      rmed=rs(n2)
      if (mod(n,2) .eq. 0) rmed=(rs(n2)+rs(n2-1))/2.0
      nup=0
      ndown=0
      nrun=1
      rp=a(1)-rmed
      do 10 i=2,n
        if (rp*(a(i)-rmed) .gt. 0.0) go to 11
        nrun=nrun+1
11      if (rp .gt. 0.0) nup=nup+1
        if (rp .le. 0.0) ndown=ndown+1
        rp=a(i)-rmed
10    continue
      if (mod(n,2) .eq. 1) nrun=nrun-1
      if (rp .gt. 0.0) nup=nup+1
      if (rp .le. 0.0) ndown=ndown+1
      return
      end
      subroutine datprt
c#    MMC routine 080 lstmod: 01/13/93
c*****Prints the date and the time
      real*8 dat,tim
      integer itim(14)
      call system('hostname')
      call system('date')
C@VX      character*8 timea
C@VX      character*9 datea
      write (6,1001)
C@BM      call datimx(itim)
C@BM      write (6,1000) itim(7),itim(6),itim(8),itim(5),itim(4),itim(3)
C@VX      call time(timea)
C@VX      call date(datea)
C@VX      write (6,1002) datea,timea
C@RY      crtime=clock()
C@RY      crdate=date()
C@RY      write (6,1002) crdate,crtime
      return
1000  format(' Date:',i4,'/',i2,'/',i4,6x,'time:',i4,':',i2,':',i2)
1001  format(' Common blocks were last modified on 08/04/98 and the',
     -  ' last program modification was done on 03/03/03',/)
1002  format(' Date: ',a9,'      Time: ',a8)
      end
