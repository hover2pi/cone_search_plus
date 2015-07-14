
c--------------------------------------------------------
c
c query_loc_2MASS.F
c
c This routine will allow super-fast searching
c of a condensed 2MASS data base.  The condensed
c data base contains RA, DEC, MJ, MH, MK and the
c quality flags for each band.  It has 13 bytes
c for each star.  With 475M stars, the total size
c of the database is 5.7G.  It can be compressed to
c fit on a single DVD for storage or export.
c
c This routine will hardcode the routine to look
c for the data base on central 
c
c Written by Jay Anderson, March 2011
c
c The original data comes from the websites associated
c with: 
c
c     SKRUTSKIE et al 2006 AJ 131 1163
c
c
c--------------------------------------------------------


c---------------------------------
c
c Set to TRUE if it's an INTEL mac,
c or something similar.  This is
c false for old MACs and SOLARIS.
c It has to do with the order that
c the machine stores integers (big
c endian/little endian).
c
#define _LINUX_ .true.
#define _DIAG_  .false.

c--------------------------------------
c
c
c
      program query_loc_2MASS
      implicit none

      integer i0
      character*80 stringg      
      character*80 STRING

      integer nstar

      real*8 ra      
      real*8 de
      real jmag, jerr
      real hmag, herr
      real kmag, kerr

c---------------------------------------
c
c variables to hold the search-qualifying 
c stars before we decide whether to output 
c them
c

      integer     U
      integer     Us

      real*8 string2ra
      real*8 string2de
      character*11 ra2sexig
      character*12 de2sexig
      character*80 RA_STRING
      character*80 DE_STRING
      character*80 STRING_IN


c-----------------------------------
c
c read in one record at a time; this
c contains the current 277-star record
c
      byte buff13_277(13,277)

      integer N

c-------------------------------------
c
c this map tells you what the first record
c is for any point in the sky
c
      integer NMAP_IJ(3601,1801)


c-------------------------------------
c
c these are the parameters that define
c the stars that quality for the search
c
c     RELATIVE-TYPE SEARCHES
c
      real*8 RCEN,  DCEN       ! (deg) the pivot RA/DEC (for radial searches)
      real   DRMIN, DRMAX      ! (deg) for radial searches or directed lists
      real   DXMAX, DYMAX      ! (deg) for localized box searches

      real   JCEN, HCEN, KCEN  ! the pivot JHK mags
      real   DJMIN, DJMAX      ! J mag lower and upper bounds (relative to list input)
      real   DHMIN, DHMAX      ! H mag ...
      real   DKMIN, DKMAX      ! K mag ...

c
c     ABSOLUTE-TYPE SEARCHES
c
      real*8 RMIN,  RMAX       ! (deg) min/max RA     (for box
      real*8 DMIN,  DMAX       ! (deg) min/max DEC     searches)
      real*8 EBMIN, EBMAX      ! (deg) min/max ecliptic longitude   (for box 
      real*8 ELMIN, ELMAX      ! (deg) min/max ecliptic latitude     searches)
      real*8 GBMIN, GBMAX      ! (deg) min/max galactic longitude   (for box
      real*8 GLMIN, GLMAX      ! (deg) min/max galactic latitude     searches)

      real   JMIN, JMAX        ! absolute J mag lower/upper bounds
      real   HMIN, HMAX        ! absolute H mag lower/upper bounds
      real   KMIN, KMAX        ! absolute K mag lower/upper bounds


      real*8      dd_lcky
      real*8      ra_lcky, de_lcky
      real*8      mj_lcky, mh_lcky, mk_lcky
      character*3 aa_lcky
      real*8      dx_lcky, dy_lcky


      logical LCKY
      logical LIST

c---------------------------------------------------------------

      integer IIMINu,IIMAXu
      integer IIMIN, IIMAX
      integer JJMIN, JJMAX

      integer i, j


      real*8 rminij, rmaxij
      real*8 dminij, dmaxij

      byte usesky_x01( 360, 180) 
      byte usesky_c01( 360, 180) 
      byte usesky_x10(3600,1800) 
      byte usesky_c10(3600,1800) 
      byte use
      integer NUSE 

      integer L, Ls
      integer jl(999) 
      integer i1l(999) 
      integer i2l(999) 

      integer LL, LLs
      integer dn_ll(9999)
      integer jj_ll(9999)
      integer ii1_ll(9999)
      integer ii2_ll(9999)
      integer nn1_ll(9999)
      integer nn2_ll(9999)

      integer rec0, rec1, rec2, recu
      integer NN, nn1, nn2

      real*8 dx, dy, dd, rd2x, rd2y
      integer ii, ii1, ii2
      integer jj, jj1, jj2
      integer ii1u, ii2u
      integer iiu,  jju

      integer iii, jjj, iiiu, jjju
c
c---------------------------------------------------------------
c
      integer*4   ira, ide
      integer*2   imj, imh, imk
      character*3 aaa

      integer NARG, NARGs
      character*1 UorV

      RCEN = -999.9
      DCEN = -999.9

      RMIN = -999.9
      RMAX = 9999.9
      DMIN = -999.9
      DMAX = 9999.9

      EBMIN = -999.9
      EBMAX = 9999.9
      ELMIN = -999.9
      ELMAX = 9999.9
       
      GBMIN = -999.9
      GBMAX = 9999.9
      GLMIN = -999.9
      GLMAX = 9999.9
      
      JMIN =  -9.9      
      JMAX =  99.9      
      HMIN =  -9.9      
      HMAX =  99.9      
      KMIN =  -9.9      
      KMAX =  99.9      

      JCEN = 0.00
      HCEN = 0.00
      KCEN = 0.00
 
      DRMIN = -999.99
      DRMAX =  999.99
      DJMIN = -999.99
      DJMAX =  999.99
      DHMIN = -999.99
      DHMAX =  999.99
      DKMIN = -999.99
      DKMAX =  999.99
      DXMAX =  999.99
      DYMAX =  999.99
      
      LCKY = .false.
      LIST = .false.
      UorV = 'U'

      NARGs = iargc()
      if (NARGs.eq.0) then
         print*,' '
         print*,'ROUTINE NEEDS COMMAND-LINE ARGUMENTS: '
         print*,' '
         print*,' query_2mass.e RA_MIN=5.0 ... '
         print*,' '
         print*,'    Specify any, or all of:    '
         print*,' '
         print*,'    RA_MIN: minimum RA (in degrees), '
         print*,'    RA_MAX: maximum RA (in degrees)  '
         print*,'    DE_MIN: minimum Dec (in degrees) '
         print*,'    DE_MAX: maximum Dec (in degrees) '
         print*,' '
         print*,'    RA_CEN: center RA point '
         print*,'    DE_CEN: center Dec point '
         print*,' '
         print*,'    DR_MIN: minimum offset from RA/DEC pt (degs)'
         print*,'    DR_MAX: maximum offset from RA/DEC pt (degs)'
         print*,'    DJ_MIN: minimum mag difference J from inp list'
         print*,'    DJ_MAX: maximum mag difference J from inp list'
         print*,'    DH_MIN: minimum mag difference H from inp list'
         print*,'    DH_MAX: maximum mag difference H from inp list'
         print*,'    DK_MIN: minimum mag difference K from inp list'
         print*,'    DK_MAX: maximum mag difference K from inp list'
         print*,' '
         print*,'    M[J,H,K]_MIN: minimum [j,h,k] mag to qualify '
         print*,'    M[J,H,K]_MAX: maximum [j,h,k] mag  to qualify'
         print*,' '
         print*,'    GB_MIN: minimum Galactic long (in degrees) '
         print*,'    GL_MIN: minimum Galactic long (in degrees) '
         print*,'    EB_MIN: minimum ecliptic long (in degrees) '
         print*,'    EL_MIN: minimum ecliptic long (in degrees) '
         print*,' '
         print*,'    LUCKY+ : returns the closest qualifying '
         print*,'             point to RA_CEN, DE_CEN        '
         print*,' '
         print*,'    RDLIST+: takes a list of RAs and DECs for'
         print*,'             input (from stdin) and searches '
         print*,'             about each of them.'
         stop
         endif 

      write(*,'(''#NARGs: '',i2)') NARGs 
      do NARG = 1, NARGs
         call getarg(NARG,STRING)
         if (STRING(1:7).eq.'RA_CEN=') then
             RCEN = string2ra(string(8:80))
            !read(STRING(8:80),*) RCEN
             write(*,'(''#NARG'',i2.2,'': '',a20,3x,a12,f12.6)') 
     .             NARG,STRING(1:20),'--->   RCEN: ',RCEN
             goto 1
             endif
         if (STRING(1:7).eq.'DE_CEN=') then
             DCEN = string2de(string(8:80))
            !read(STRING(8:80),*) DCEN
             write(*,'(''#NARG'',i2.2,'': '',a20,3x,a12,f12.6)') 
     .             NARG,STRING(1:20),'--->   DCEN: ',DCEN
             goto 1
             endif
         if (STRING(1:7).eq.'RA_MIN=') then
             RMIN = string2ra(string(8:80))
            !read(STRING(8:80),*) RMIN
             write(*,'(''#NARG'',i2.2,'': '',a20,3x,a12,f12.6)') 
     .             NARG,STRING(1:20),'--->   RMIN: ',RMIN
             goto 1
             endif
         if (STRING(1:7).eq.'RA_MAX=') then
             RMAX = string2ra(string(8:80))
            !read(STRING(8:80),*) RMAX
             write(*,'(''#NARG'',i2.2,'': '',a20,3x,a12,f12.6)') 
     .             NARG,STRING(1:20),'--->   RMAX: ',RMAX
             goto 1
             endif
         if (STRING(1:7).eq.'DE_MIN=') then
             DMIN = string2de(string(8:80))
            !read(STRING(8:80),*) DMIN
             write(*,'(''#NARG'',i2.2,'': '',a20,3x,a12,f12.6)') 
     .             NARG,STRING(1:20),'--->   DMIN: ',DMIN
             goto 1
             endif
         if (STRING(1:7).eq.'DE_MAX=') then
             DMAX = string2de(string(8:80))
            !read(STRING(8:80),*) DMAX
             write(*,'(''#NARG'',i2.2,'': '',a20,3x,a12,f12.6)') 
     .             NARG,STRING(1:20),'--->   DMAX: ',DMAX
             goto 1
             endif
         if (STRING(1:7).eq.'EB_MIN=') then
             read(STRING(8:80),*) EBMIN
             goto 1
             endif
         if (STRING(1:7).eq.'EB_MAX=') then
             read(STRING(8:80),*) EBMAX
             goto 1
             endif
         if (STRING(1:7).eq.'EL_MIN=') then
             read(STRING(8:80),*) ELMIN
             goto 1
             endif
         if (STRING(1:7).eq.'EL_MAX=') then
             read(STRING(8:80),*) ELMAX
             goto 1
             endif
         if (STRING(1:7).eq.'GB_MIN=') then
             read(STRING(8:80),*) GBMIN
             goto 1
             endif
         if (STRING(1:7).eq.'GB_MAX=') then
             read(STRING(8:80),*) GBMAX
             goto 1
             endif
         if (STRING(1:7).eq.'GL_MIN=') then
             read(STRING(8:80),*) GLMIN
             goto 1
             endif
         if (STRING(1:7).eq.'GL_MAX=') then
             read(STRING(8:80),*) GLMAX
             goto 1
             endif
         if (STRING(1:7).eq.'MJ_MIN=') then
             read(STRING(8:80),*) JMIN
             goto 1
             endif
         if (STRING(1:7).eq.'MJ_MAX=') then
             read(STRING(8:80),*) JMAX
             goto 1
             endif
         if (STRING(1:7).eq.'MH_MIN=') then
             read(STRING(8:80),*) HMIN
             goto 1
             endif
         if (STRING(1:7).eq.'MH_MAX=') then
             read(STRING(8:80),*) HMAX
             goto 1
             endif
         if (STRING(1:7).eq.'MK_MIN=') then
             read(STRING(8:80),*) KMIN
             goto 1
             endif
         if (STRING(1:7).eq.'MK_MAX=') then
             read(STRING(8:80),*) KMAX
             goto 1
             endif
         if (STRING(1:7).eq.'DR_MIN=') then
             read(STRING(8:80),*) DRMIN
             goto 1
             endif
         if (STRING(1:7).eq.'DR_MAX=') then
             read(STRING(8:80),*) DRMAX
             goto 1
             endif
         if (STRING(1:7).eq.'DJ_MIN=') then
             read(STRING(8:80),*) DJMIN
             goto 1
             endif
         if (STRING(1:7).eq.'DJ_MAX=') then
             read(STRING(8:80),*) DJMAX
             goto 1
             endif
         if (STRING(1:7).eq.'DH_MIN=') then
             read(STRING(8:80),*) DHMIN
             goto 1
             endif
         if (STRING(1:7).eq.'DH_MAX=') then
             read(STRING(8:80),*) DHMAX
             goto 1
             endif
         if (STRING(1:7).eq.'DK_MIN=') then
             read(STRING(8:80),*) DKMIN
             goto 1
             endif
         if (STRING(1:7).eq.'DK_MAX=') then
             read(STRING(8:80),*) DKMAX
             goto 1
             endif
         if (STRING(1:7).eq.'DX_MAX=') then
             read(STRING(8:80),*) DXMAX
             goto 1
             endif
         if (STRING(1:7).eq.'DY_MAX=') then
             read(STRING(8:80),*) DYMAX
             goto 1
             endif
         if (STRING(1:6).eq.'+LUCKY'.or.
     .       STRING(1:6).eq.'LUCKY+') then
             LCKY = .true.
             goto 1
             endif
         if (STRING(1:8).eq.'+RDLIST'.or.
     .       STRING(1:8).eq.'RDLIST+') then
             LIST = .true.
             UorV = 'V'
             goto 1
             endif

         print*,'ERROR:  COULD NOT PARSE COMMAND LINE.  '
         print*,'UNRECOGNIZED FIELD SPECIFIER:'
         print*,'    STRING(1:7): ',STRING(1:7)
         print*,'    STRING(8:+): ',STRING(8:40)
         stop
    1    continue
         enddo

       if (LCKY.and.RCEN.lt.0.00.and.(.not.LIST)) then
          print*,' '
          print*,'IF YOU WANT TO DO "I FEEL LUCKY" CALCULATIONS,'
          print*,'YOU NEED TO HAVE A POINT (OR A LIST OF POINTS)'
          print*,'ABOUT WHICH THE SEARCH IS BEING CONDUCTED.'
          print*,' '
          endif


c---------------------------------------------------------
c
c this will read in the map to the database
c
      write(*,'(''# OPEN BASE.DATA for map...'')')
      open(12,file='BASE.DATA',status='old',
     .        recl= 25941604,
     .        form='UNFORMATTED',
     .        access='DIRECT')
      read(12,rec=1) NMAP_IJ
      close(12)


c---------------------------------------------------------
c
c now, open the file for actual use as the dB
c
      rec0  = 3601*1801*4 / 13 
      write(*,'(''# OPEN BASE.DATA for work...'')')
      open(12,file='BASE.DATA',status='old',
     .        recl=3601,form='UNFORMATTED',
     .        access='DIRECT')


 137  continue

      STRING_IN = ''
      if (LIST) then
         read(*,'(a80)',end=777) STRING_IN
         if (STRING_IN(1:1).eq.'#') goto 137
         read(STRING_IN,*) RA_STRING, DE_STRING, JCEN, HCEN, KCEN
         RCEN = string2ra(RA_STRING)
         DCEN = string2de(DE_STRING)
         endif

      do i = 001, 360
      do j = 001, 180
         usesky_x01(i,j) = -100
         enddo
         enddo

      i = 1 + (RCEN)      
      j = 1 + (DCEN + 90) 
      if (i.ge.001.and.i.le.360.and.
     .    j.ge.001.and.j.le.180) usesky_x01(i,j) = 1


      do i = 001, 360
      do j = 001, 180
         rminij = (i-1)*1.0 
         rmaxij = (i  )*1.0
         dminij = (j-1)*1.0 - 90.00
         dmaxij = (j  )*1.0 - 90.00
         usesky_x01(i,j) = 
     .       max(usesky_x01(i,j),
     .           use(rminij,dminij,jmin,hmin,kmin,
     .               rcen,dcen,jcen,hcen,kcen,
     .                    drmin,drmax,
     .                    djmin,djmax,dhmin,dhmax,dkmin,dkmax,
     .                    dxmax,dymax,
     .               rmin,rmax,dmin,dmax,
     .               ebmin,ebmax,elmin,elmax,
     .               gbmin,gbmax,glmin,glmax,
     .               jmin,jmax,hmin,hmax,kmin,kmax),
     .           use(rmaxij,dminij,jmin,hmin,kmin,
     .               rcen,dcen,jcen,hcen,kcen,
     .                    drmin,drmax,
     .                    djmin,djmax,dhmin,dhmax,dkmin,dkmax,
     .                    dxmax,dymax,
     .               rmin,rmax,dmin,dmax,
     .               ebmin,ebmax,elmin,elmax,
     .               gbmin,gbmax,glmin,glmax,
     .               jmin,jmax,hmin,hmax,kmin,kmax),
     .           use(rminij,dmaxij,jmin,hmin,kmin,
     .               rcen,dcen,jcen,hcen,kcen,
     .                    drmin,drmax,
     .                    djmin,djmax,dhmin,dhmax,dkmin,dkmax,
     .                    dxmax,dymax,
     .               rmin,rmax,dmin,dmax,
     .               ebmin,ebmax,elmin,elmax,
     .               gbmin,gbmax,glmin,glmax,
     .               jmin,jmax,hmin,hmax,kmin,kmax),
     .           use(rmaxij,dmaxij,jmin,hmin,kmin,
     .               rcen,dcen,jcen,hcen,kcen,
     .                    drmin,drmax,
     .                    djmin,djmax,dhmin,dhmax,dkmin,dkmax,
     .                    dxmax,dymax,
     .               rmin,rmax,dmin,dmax,
     .               ebmin,ebmax,elmin,elmax,
     .               gbmin,gbmax,glmin,glmax,
     .               jmin,jmax,hmin,hmax,kmin,kmax))
         enddo
         enddo


      do i = 001, 360
      do j = 001, 180
         usesky_c01(i,j) = 0
         do ii = i-1, i+1
         do jj = j-1, j+1
            iiu = ii
            if (iiu.eq.361) iiu = 001
            if (iiu.eq.000) iiu = 360
            jju = jj
            if (jju.eq.181) jju = 180
            if (jju.eq.000) jju = 001
            if (usesky_x01(iiu,jju).eq.1) usesky_c01(i,j) = 1
            enddo 
            enddo 
         enddo 
         enddo 

      if (_DIAG_) then
         write(*,'(''# DIAG output usesky_[x,c]01.fits'')')
         call writfits_b1('usesky_x01.fits',usesky_x01,360,180)
         call writfits_b1('usesky_c01.fits',usesky_c01,360,180)
         endif

      do i = 0001, 3600
      do j = 0001, 1800
         usesky_x10(i,j) = 0
         usesky_c10(i,j) = 0
         enddo
         enddo


      i = 1 + (RCEN   )*10
      j = 1 + (DCEN+90)*10
      if (i.ge.001.and.i.le.3600.and.
     .    j.ge.001.and.j.le.1800) usesky_x10(i,j) = 1

      if (_DIAG_) then
         print*,'  jmin: ',jmin,hmin,kmin
         print*,'  rcen: ',rcen,dcen
         print*,'  jcen: ',jcen,hcen,kcen
         print*,' drmin: ',drmin,drmax
         print*,' djmin: ',djmin,djmax
         print*,' dhmin: ',dhmin,dhmax
         print*,' dkmin: ',dkmin,dkmax
         print*,' dxmax: ',dxmax,dymax
         print*,'  rmin: ',rmin,rmax
         print*,'  dmin: ',dmin,dmax
         print*,'  jmin: ',jmin,jmax
         print*,'  hmin: ',hmin,hmax
         print*,'  kmin: ',kmin,kmax
         endif

      JJMIN =  9999
      JJMAX = -9999
      IIMIN =  9999
      IIMAX = -9999
      do i = 001, 360
      do j = 001, 180
         ii1 = 1 + (i-1)*10 -1
         ii2 = 1 + (i  )*10 +1
         jj1 = 1 + (j-1)*10 -1
         jj2 = 1 + (j  )*10 +1
         if (usesky_c01(i,j).eq.1) then
            do jj = jj1, jj2
               jju = jj
               if (jju.lt.0001) jju = 0001
               if (jju.gt.1800) jju = 1800
               do ii = ii1, ii2
                  iiu = ii
                  if (iiu.lt.0001) iiu = iiu+3600
                  if (iiu.gt.3600) iiu = iiu-3600
                  rminij = (iiu-1)*0.1
                  rmaxij = (iiu  )*0.1
                  dminij = (jju-1)*0.1 - 90.00
                  dmaxij = (jju  )*0.1 - 90.00
                  usesky_x10(iiu,jju) = 
     .                max(usesky_x10(iiu,jju),
     .                    use(rminij,dminij,jmin,hmin,kmin,
     .                        rcen,dcen,jcen,hcen,kcen,
     .                             drmin,drmax,
     .                             djmin,djmax,dhmin,dhmax,dkmin,dkmax,
     .                             dxmax,dymax,
     .                        rmin,rmax,dmin,dmax,
     .                        ebmin,ebmax,elmin,elmax,
     .                        gbmin,gbmax,glmin,glmax,
     .                        jmin,jmax,hmin,hmax,kmin,kmax),
     .                    use(rmaxij,dminij,jmin,hmin,kmin,
     .                        rcen,dcen,jcen,hcen,kcen,
     .                             drmin,drmax,
     .                             djmin,djmax,dhmin,dhmax,dkmin,dkmax,
     .                             dxmax,dymax,
     .                        rmin,rmax,dmin,dmax,
     .                        ebmin,ebmax,elmin,elmax,
     .                        gbmin,gbmax,glmin,glmax,
     .                        jmin,jmax,hmin,hmax,kmin,kmax),
     .                    use(rminij,dmaxij,jmin,hmin,kmin,
     .                        rcen,dcen,jcen,hcen,kcen,
     .                             drmin,drmax,
     .                             djmin,djmax,dhmin,dhmax,dkmin,dkmax,
     .                             dxmax,dymax,
     .                        rmin,rmax,dmin,dmax,
     .                        ebmin,ebmax,elmin,elmax,
     .                        gbmin,gbmax,glmin,glmax,
     .                        jmin,jmax,hmin,hmax,kmin,kmax),
     .                    use(rmaxij,dmaxij,jmin,hmin,kmin,
     .                        rcen,dcen,jcen,hcen,kcen,
     .                             drmin,drmax,
     .                             djmin,djmax,dhmin,dhmax,dkmin,dkmax,
     .                             dxmax,dymax,
     .                        rmin,rmax,dmin,dmax,
     .                        ebmin,ebmax,elmin,elmax,
     .                        gbmin,gbmax,glmin,glmax,
     .                        jmin,jmax,hmin,hmax,kmin,kmax))
                  if (usesky_x10(iiu,jju).eq.1) then
                     do iii = iiu-1, iiu+1
                     do jjj = jju-1, jju+1
                        iiiu = iii
                        jjju = jjj
                        if (iiiu.lt.0001) iiiu = iiiu+3600
                        if (iiiu.gt.3600) iiiu = iiiu-3600
                        if (jjju.lt.0001) jjju = 0001
                        if (jjju.gt.1800) jjju = 1800
                        usesky_c10(iiiu,jjju) = 1
                        if (iiiu.lt.IIMIN) IIMIN = iiiu
                        if (iiiu.gt.IIMAX) IIMAX = iiiu
                        if (jjju.lt.JJMIN) JJMIN = jjju
                        if (jjju.gt.JJMAX) JJMAX = jjju
                        enddo
                        enddo
                     endif
                  enddo  
                enddo
             endif
          enddo
          enddo

      if (_DIAG_) then 
         print*,' '
         print*,'REGION OF SKY TO SEARCH: '
         print*,' '
         print*,'  IIMIN/MAX: ',IIMIN,IIMAX
         print*,'  JJMIN/MAX: ',JJMIN,JJMAX
         print*,' '
         endif

      if (_DIAG_) then 
         write(*,'(''# DIAG output usesky_[x,c]10.fits'')')
         call writfits_b1('usesky_x10.fits',usesky_x10,3600,1800)
         call writfits_b1('usesky_c10.fits',usesky_c10,3600,1800)
         endif


      open(15,file='query_2mass.out',status='unknown')

      dd_lcky = 999.999
      ra_lcky = 0.00
      de_lcky = 0.00
      mj_lcky = 99.99
      mh_lcky = 99.99
      mk_lcky = 99.99
      aa_lcky = '???'
      dx_lcky = 99.99
      dy_lcky = 99.99
      Us = 0
      do jj = JJMIN, JJMAX
         IIMINu = IIMAX
         IIMAXu = IIMIN
         do ii = IIMIN, IIMAX
            if (usesky_c10(ii,jj).eq.1) then
               if (ii.lt.IIMINu) IIMINu = ii
               if (ii.gt.IIMAXu) IIMAXu = ii
               endif
            enddo
         nn1 = max(NMAP_IJ(IIMINu  ,jj),rec0+1)
         nn2 = max(NMAP_IJ(IIMAXu+1,jj),rec0+1)
         rec1 = nn1/277 + 1
         rec2 = nn2/277 + 1
         do N = rec1, rec2
            read(12,rec=N) buff13_277
            do NN = 1, 277
               call buff2rdm(buff13_277(1,NN),ira,ide,
     .                                        imj,imh,imk,aaa)
               ra = ira/60./60./100.
               de = -90 + ide/60./60./100.
               jmag = imj/100.
               hmag = imh/100.
               kmag = imk/100.
               U = use(ra,de,jmag,hmag,kmag,
     .                 rcen,dcen,jcen,hcen,kcen,
     .                      drmin,drmax,
     .                      djmin,djmax,dhmin,dhmax,dkmin,dkmax,
     .                      dxmax,dymax,
     .                 rmin,rmax,dmin,dmax,
     .                 ebmin,ebmax,elmin,elmax,
     .                 gbmin,gbmax,glmin,glmax,
     .                 jmin,jmax,hmin,hmax,kmin,kmax)
               if (U.eq.1) then
                  dx = 0.
                  dy = 0.
                  if (RCEN.gt.00) then
                     dx = rd2x(ra,de,rcen,dcen)
                     dy = rd2y(ra,de,rcen,dcen)
                     dd = sqrt(dx**2+dy**2)
                     if (dd.lt.dd_lcky) then
                        dd_lcky = dd
                        ra_lcky = ra
                        de_lcky = de
                        mj_lcky = jmag
                        mh_lcky = hmag
                        mk_lcky = kmag
                        aa_lcky = aaa
                        dx_lcky = dx
                        dy_lcky = dy
                        endif
                     endif
                  Us = Us + 1
                  if (.not.LCKY) then
                     if (RCEN.ge.0) then
                         write(15,159) ra,de,
     .                                 jmag,hmag,kmag,
     .                                 aaa,UorV,Us,
     .                                 dx,dy,STRING_IN
                         write( *,159) ra,de,
     .                                jmag,hmag,kmag,
     .                                aaa,UorV,Us,
     .                                dx,dy,STRING_IN
                         endif
                     if (RCEN.lt.0) then
                         write(15,159) ra,de,
     .                                 jmag,hmag,kmag,
     .                                 aaa,UorV,Us
                         write( *,159) ra,de,
     .                                 jmag,hmag,kmag,
     .                                 aaa,UorV,Us
                         endif
                     endif
                  endif
               enddo!NN
            enddo!N
         enddo!jj


  159 format(2x,f9.5,1x,f9.5,1x,3f8.3,1x,a3,
     .       2x,a1,i8.8,
     .       2x,2(f9.5,1x),
     .       2x,80a)


      if (LCKY) then
         print*,' '
         print*,'"LUCKY" RESULT: '
         print*,' '
         write(*,159) ra_lcky,de_lcky,
     .                mj_lcky,mh_lcky,mk_lcky,
     .                aa_lcky,UorV,Us,
     .                dx_lcky,dy_lcky,
     .                STRING_IN
         print*,' '
         endif


      if (LIST) goto 137
      close(12)
      close(15)

      stop

 777  continue
      print*,' '
      print*,'END OF INPUT FILE...'
      print*,' '

      stop
      end




c--------------------------------------------------------------------
c
c
      byte function use(r,d,mj,mh,mk,
     .                  rcen,dcen,jcen,hcen,kcen,
     .                       drmin,drmax,
     .                       djmin,djmax,dhmin,dhmax,dkmin,dkmax,
     .                       dxmax,dymax,
     .                  rmin,rmax,dmin,dmax,
     .                  ebmin,ebmax,elmin,elmax,
     .                  gbmin,gbmax,glmin,glmax,
     .                  jmin,jmax,hmin,hmax,kmin,kmax)
      implicit none

      real*8 r,d
      real   mj,mh,mk

      real*8 rcen, dcen
      real   jcen, hcen, kcen
      real   drmin, drmax
      real   djmin, djmax, dhmin, dhmax, dkmin, dkmax
      real   dxmax, dymax
      real*8  rmin,  rmax,  dmin,  dmax
      real*8 ebmin, ebmax, elmin, elmax
      real*8 gbmin, gbmax, glmin, glmax
      real   jmin, jmax
      real   hmin, hmax
      real   kmin, kmax

      real*8 gl, gb
      real*8 el, eb

      real*8 rd2x, rd2y
      real*8 dx, dy, dd

      if (rmin.gt.0) then
         use = -1
         if (r.lt.rmin) return
         use = -2
         if (r.gt.rmax) return
         use = -3
         if (d.lt.dmin) return
         use = -4
         if (d.gt.dmax) return
         endif
 
      if (ebmin.gt.-99.9) then
         call rd2lb_E(r,d,el,eb)
         use = -5
         if (el.lt.elmin) return
         use = -6
         if (el.gt.elmax) return
         use = -7
         if (eb.lt.ebmin) return
         use = -8
         if (eb.gt.ebmax) return
         endif
      
      if (gbmin.gt.-99.9) then
         call rd2lb_G(r,d,gl,gb)
         use = -9
         if (gl.lt.glmin) return
         use = -10
         if (gl.gt.glmax) return
         use = -11
         if (gb.lt.gbmin) return
         use = -12
         if (gb.gt.gbmax) return
         endif
    
      if (jmin.gt.-9.99) then
         use = -13
         if (mj.lt.jmin) return
         use = -14
         if (mj.gt.jmax) return
         endif

      if (hmin.gt.-9.99) then
         use = -15
         if (mh.lt.hmin) return
         use = -16
         if (mh.gt.hmax) return
         endif

      if (kmin.gt.-9.99) then
         use = -17
         if (mk.lt.kmin) return
         use = -18
         if (mk.gt.kmax) return
         endif

      if (rcen.ge.0.00) then
         dx = rd2x(r,d,rcen,dcen)
         dy = rd2y(r,d,rcen,dcen)
         dd = sqrt(dx**2+dy**2)
         use = -19
         if (dd.lt.drmin) return
         use = -20
         if (dd.gt.drmax) return
         use = -21
         if (mj-jcen.gt.djmax) return
         use = -22
         if (mj-jcen.lt.djmin) return
         use = -23
         if (mh-hcen.gt.dhmax) return
         use = -24
         if (mh-hcen.lt.dhmin) return
         use = -25
         if (mk-kcen.gt.dkmax) return
         use = -26
         if (mk-kcen.lt.dkmin) return
         use = -27
         if (abs(dx).gt.dxmax) return
         use = -28
         if (abs(dy).gt.dymax) return
         endif

      use = 1

      return
      end



c-----------------------------------------------
c 
c This routine will take an (ra,dec) and a 
c tangent-point (ra0,dec0) and will compute
c the x coordinate in the tangent plane.  The 
c tangent-plane x-axis is directed along -RA.
c The units for x will be degrees, just lilke
c those for RA and DEC.
c
      real*8 function rd2x(r,d,r0,d0)
      implicit none

      real*8 r, d
      real*8 r0,d0

      real*8 cosra, sinra
      real*8 cosde, sinde
      real*8 cosd0, sind0
      real*8 rrrr
      real*8 xrad
      real   x, y, z
      real   xx,yy,zz

      cosra = cos((r-r0)*3.141592654d0/180.0d0)
      sinra = sin((r-r0)*3.141592654d0/180.0d0)

      cosde = cos(d *3.141592654d0/180.0d0)
      sinde = sin(d *3.141592654d0/180.0d0)
      cosd0 = cos(d0*3.141592654d0/180.0d0)
      sind0 = sin(d0*3.141592654d0/180.0d0)
  
      rrrr = sind0*sinde + cosd0*cosde*cosra


      xrad = cosde*sinra/rrrr
      rd2x = xrad*180.0d0/3.141592654d0
  
      x  = cosde*cos(r *3.14159/180)
      y  = cosde*sin(r *3.14159/180)
      z  = sinde
      xx = cosd0*cos(r0*3.14159/180)
      yy = cosd0*sin(r0*3.14159/180)
      zz = sind0

      if (x*xx + y*yy + z*zz.lt.0) rd2x = 90


      return
      end


c-----------------------------------------------
c 
c This routine will take an (ra,dec) and a 
c tangent-point (ra0,dec0) and will compute
c the y coordinate in the tangent plane.  The 
c tangent-plane y-axis is directed along +DEC.
c The units for y will be degrees, just lilke
c those for RA and DEC.
c
      real*8 function rd2y(r,d,r0,d0)
      implicit none

      real*8 r, d
      real*8 r0,d0

      real*8 cosra, sinra
      real*8 cosde, sinde
      real*8 cosd0, sind0
      real*8 rrrr
      real*8 yrad
      real   x, y, z
      real   xx,yy,zz

      cosra = cos((r-r0)*3.141592654d0/180.0d0)
      sinra = sin((r-r0)*3.141592654d0/180.0d0)
      cosde = cos(d *3.141592654d0/180.0d0)
      sinde = sin(d *3.141592654d0/180.0d0)
      cosd0 = cos(d0*3.141592654d0/180.0d0)
      sind0 = sin(d0*3.141592654d0/180.0d0)

      rrrr = sind0*sinde + cosd0*cosde*cosra

      yrad = (cosd0*sinde-sind0*cosde*cosra)/rrrr
      rd2y = yrad*180.0d0/3.141592654d0

      x  = cosde*cos(r *3.14159/180)
      y  = cosde*sin(r *3.14159/180)
      z  = sinde
      xx = cosd0*cos(r0*3.14159/180)
      yy = cosd0*sin(r0*3.14159/180)
      zz = sind0

      if (x*xx + y*yy + z*zz.lt.0) rd2y = 90

      return
      end


c-----------------------------------------------
c 
c This routine will take a point (x,y) in the
c tangent plane, and the coordinate of the
c tangent-point (ra0,dec0), and will compute
c the ra coordinate that corresponds to (x,y).
c
      real*8 function xy2r(x,y,r0,d0)
      implicit none

      real*8 x, y
      real*8 r0,d0

      real*8 cosde, sinde
      real*8 cosd0, sind0
      real*8 tandr, dr
      real*8 xrad, yrad

      xrad = x/180.0d0*3.141592654d0
      yrad = y/180.0d0*3.141592654d0

      cosd0 = cos(d0*3.141592654d0/180.0d0)
      sind0 = sin(d0*3.141592654d0/180.0d0)

      tandr = xrad/(cosd0-yrad*sind0)
      dr    = atan(tandr)
      xy2r  = r0 + dr*180.0d0/3.141592654d0

      return
      end




c-----------------------------------------------
c 
c This routine will take a point (x,y) in the
c tangent plane, and the coordinate of the
c tangent-point (ra0,dec0), and will compute
c the ra coordinate that corresponds to (x,y).
c
      real*8 function xy2d(x,y,r0,d0)
      implicit none

      real*8 x, y
      real*8 r0,d0

      real*8 cosd0, sind0
      real*8 tandr, dr
      real*8 cosdr
      real*8 tande
      real*8 xrad, yrad

      xrad = x/180.0d0*3.141592654d0
      yrad = y/180.0d0*3.141592654d0

      cosd0 = cos(d0*3.141592654d0/180.0d0)
      sind0 = sin(d0*3.141592654d0/180.0d0)

      tandr = xrad/(cosd0-yrad*sind0)
      dr    = atan(tandr)

      cosdr = cos(dr)

      tande = cosdr*(sind0+yrad*cosd0)/
     .              (cosd0-yrad*sind0)

      print*,'---> tande: ',tande 
      xy2d = atan(tande)*180.0d0/3.141592654d0

      return
      end


 
 
c--------------------------------------
c
c ra,dec to Galactic longitude, latitude
c
      subroutine rd2lb_G(ra,dec,l,b)
      implicit none
 
      real*8 ra, dec
      real*8 l, b
 
      b = 180/3.14159*asin(cos(3.14159/180*dec)
     .                    *cos(3.14159/180*27.4)
     .                    *cos(3.14159/180*(ra-192.25))
     .                    +sin(3.14159/180*dec)
     .                    *sin(3.14159/180*27.4))
 
      l = 180/3.14159*atan2(sin(3.14159/180*dec)
     .                     -sin(3.14159/180*b)
     .                     *sin(3.14159/180*27.4)
     .                     ,cos(3.14159/180*dec)
     .                     *sin(3.14159/180*(ra-192.25))
     .                     *cos(3.14159/180*27.4)) + 33.
 
      if (b.gt.180) b = b - 360 
      if (l.lt.000) l = l + 360
 
      return
      end
 
 
 
 
c--------------------------------------
c
c Galactic longitude, latitude to ra, dec
c
      subroutine lb2rd_G(l,b,ra,dec)
      implicit none
 
      real*8 l, b
      real*8 ra, dec
 
      dec = 180/3.14159*asin(cos(3.14159/180*b)
     .                      *cos(3.14159/180*27.4)
     .                      *sin(3.14159/180*(l-33.0))
     .                      +sin(3.14159/180*b)
     .                      *sin(3.14159/180*27.4))
 
      ra = 180/3.14159*atan2(cos(3.14159/180*b)
     .                      *cos(3.14159/180*(l-33))
     .                      ,sin(3.14159/180*b)
     .                      *cos(3.14159/180*27.4)
     .                      -cos(3.14159/180*b)
     .                      *sin(3.14159/180*27.4)
     .                      *sin(3.14159/180*(l-33))) + 192.25
 
      if (l.lt.0) l = l + 360
 
      return
      end
 
 
c--------------------------------------
c
c ra,dec to ecliptic longitude, latitude
c
      subroutine rd2lb_E(ra,dec,l,b)
      implicit none
 
      real*8 ra, dec
      real*8 l, b
 
      b = 180/3.14159*asin(sin(3.14159/180*dec)
     .                    *cos(3.14159/180*23.441)
     .                    -cos(3.14159/180*dec)
     .                    *sin(3.14159/180*23.441)
     .                    *sin(3.14159/180*ra))
 
      l = 180/3.14159*atan2(sin(3.14159/180*ra)
     .                     *cos(3.14159/180*23.411)
     .                     +tan(3.14159/180*dec)
     .                     *sin(3.14159/180*23.411)
     .                     ,cos(3.14159/180*ra))

      if (b.gt.180) b = b - 360 
      if (l.lt.000) l = l + 360
 
      return
      end
 
 
c--------------------------------------
c
c ra,dec to ecliptic longitude, latitude
c
      subroutine lb2rd_E(l,b,ra,dec)
      implicit none
 
      real*8 l, b
      real*8 ra, dec
 
      dec = 180/3.14159*asin(sin(3.14159/180*23.441)
     .                      *sin(3.14159/180*l)
     .                      *cos(3.14159/180*b)
     .                      +cos(3.14159/180*23.441)
     .                      *sin(3.14159/180*b))
 
      ra  = 180/3.14159*atan2(cos(3.14159/180*23.441)
     .                       *sin(3.14159/180*l)
     .                       *cos(3.14159/180*b)
     .                       -sin(3.14159/180*23.411)
     .                       *sin(3.14159/180*b)
     .                       ,cos(3.14159/180*l)
     .                       *cos(3.14159/180*b))
 
      if (ra.lt.0) ra = ra + 360
 
      return
      end


 
      subroutine buff2rdmX(buff_12,ra,de,jmag,hmag,kmag,QUAL)
      implicit none

      byte        buff_12(12)
      real*8      ra, de
      real        jmag, hmag, kmag
      character*3 QUAL

      integer     b

      integer ira, irao
      integer ide, ideo
      integer*2 imj, jqual
      integer*2 imh, hqual
      integer*2 imk, kqual

      byte        buffb_i4(4)
      integer*4   buffi_i4
      equivalence (buffb_i4,buffi_i4)

      byte        buffb_i2(2)
      integer*2   buffi_i2
      equivalence (buffb_i2,buffi_i2)

      byte   buffb_ra(3)
      byte   buffb_de(3)
      byte   buffb_mj(2)
      byte   buffb_mh(2)
      byte   buffb_mk(2)


      irao = 0
      ideo = 0

      buffb_ra(01) = buff_12(01)
      buffb_ra(02) = buff_12(02)
      buffb_ra(03) = buff_12(03)
      buffb_de(01) = buff_12(04)
      buffb_de(02) = buff_12(05)
      buffb_de(03) = buff_12(06)
      buffb_mj(01) = buff_12(07)
      buffb_mj(02) = buff_12(08)
      buffb_mh(01) = buff_12(09)
      buffb_mh(02) = buff_12(10)
      buffb_mk(01) = buff_12(11)
      buffb_mk(02) = buff_12(12)

cc---------------------------------
      buffb_i4(1) = buffb_ra(1)
      buffb_i4(2) = buffb_ra(2)
      buffb_i4(3) = buffb_ra(3)
      buffb_i4(4) = 0
      ira = buffi_i4 

cc---------------------------------
      buffb_i4(1) = buffb_de(1)
      buffb_i4(2) = buffb_de(2)
      buffb_i4(3) = buffb_de(3)
      buffb_i4(4) = 0
      ide = buffi_i4 

cc---------------------------------
      buffb_i2(1) = buffb_mj(1)
      buffb_i2(2) = buffb_mj(2)
      imj   = buffi_i2/8

      jqual = buffi_i2 - buffi_i2/8*8
      if (jqual.eq.0) QUAL(1:1) = 'X'
      if (jqual.eq.1) QUAL(1:1) = 'U'
      if (jqual.eq.2) QUAL(1:1) = 'F'
      if (jqual.eq.3) QUAL(1:1) = 'E'
      if (jqual.eq.4) QUAL(1:1) = 'A'
      if (jqual.eq.5) QUAL(1:1) = 'B'
      if (jqual.eq.6) QUAL(1:1) = 'C'
      if (jqual.eq.7) QUAL(1:1) = 'D'


cc---------------------------------
      buffb_i2(1) = buffb_mh(1)
      buffb_i2(2) = buffb_mh(2)
      imh = buffi_i2/8

      hqual = buffi_i2 - buffi_i2/8*8
      if (hqual.eq.0) QUAL(2:2) = 'X'
      if (hqual.eq.1) QUAL(2:2) = 'U'
      if (hqual.eq.2) QUAL(2:2) = 'F'
      if (hqual.eq.3) QUAL(2:2) = 'E'
      if (hqual.eq.4) QUAL(2:2) = 'A'
      if (hqual.eq.5) QUAL(2:2) = 'B'
      if (hqual.eq.6) QUAL(2:2) = 'C'
      if (hqual.eq.7) QUAL(2:2) = 'D'

cc---------------------------------
      buffb_i2(1) = buffb_mk(1)
      buffb_i2(2) = buffb_mk(2)
      imk = buffi_i2/8

      kqual = buffi_i2 - buffi_i2/8*8
      if (kqual.eq.0) QUAL(3:3) = 'X'
      if (kqual.eq.1) QUAL(3:3) = 'U'
      if (kqual.eq.2) QUAL(3:3) = 'F'
      if (kqual.eq.3) QUAL(3:3) = 'E'
      if (kqual.eq.4) QUAL(3:3) = 'A'
      if (kqual.eq.5) QUAL(3:3) = 'B'
      if (kqual.eq.6) QUAL(3:3) = 'C'
      if (kqual.eq.7) QUAL(3:3) = 'D'

      ra   = ira/10000.0
      de   = ide/10000.0 - 90.00
      jmag = imj/100.
      hmag = imh/100.
      kmag = imk/100.

      end




c---------------------------------------------------
c
      integer function ubyte(b)
      implicit none

      byte   b
      ubyte = b
      if (ubyte.lt.0) ubyte = ubyte + 256

      return
      end


c---------------------------------------------------
c
c designed to compress better...
c
c
      subroutine rdm2buff(ra,de,jmag,hmag,kmag,QUAL,buff13)
      implicit none

      real*8      ra, de
      real        jmag, hmag, kmag
      character*3 QUAL
      byte        buff13(13)

      integer     b

      integer*4 ira, ira1, ira2, iras
      integer*4 ide, ide1, ide2, ides
      integer*2 imj, jqual
      integer*2 imh, hqual
      integer*2 imk, kqual

      integer*4 irax
      integer*4 idex

      integer*2 jmags, jmagx
      integer*2 hmags, hmagx
      integer*2 kmags, kmagx

      byte b32_rai(32)
      byte b32_dei(32)
      byte b16_mji(16)
      byte b16_mhi(16)
      byte b16_mki(16)

      integer*4 ira0
      integer*4 ide0
      integer*2 ijm0
      integer*2 ihm0
      integer*2 ikm0

      ira = ((ra   )*60*60*100) + 0.5 ! max: 129,600,000 :  2**27 = 134,217,728  ! 8+8+8+3 bytes
      ide = ((de+90)*60*60*100) + 0.5 ! max:  64,800,000 :  2**26 =  67,108,864  ! 8+8+8+2 bytes


      ira1 =  ira/36000
      ira2 = (ira-ira/36000*36000)

      ide1 =  ide/36000
      ide2 = (ide-ide/36000*36000)

      iras = ira1*65536 + ira2
      ides = ide1*65536 + ide2

      jmags = abs(jmag*100) + 0.5
      if (jmags.gt. 2000) jmags =  2000
      if (jmag.lt.0) jmags = -jmags

      hmags = abs(hmag*100) + 0.5
      if (hmags.gt. 2000) hmags =  2000
      if (hmag.lt.0) hmags = -hmags

      kmags = abs(kmag*100) + 0.5
      if (kmags.gt. 2000) kmags =  2000
      if (kmag.lt.0) kmags = -kmags

      ira0 = ira
      ide0 = ide
      ijm0 = jmags
      ihm0 = hmags
      ikm0 = kmags

      jqual = -1
      if (QUAL(1:1).eq.'X') jqual = 0
      if (QUAL(1:1).eq.'U') jqual = 1
      if (QUAL(1:1).eq.'F') jqual = 2
      if (QUAL(1:1).eq.'E') jqual = 3
      if (QUAL(1:1).eq.'A') jqual = 4
      if (QUAL(1:1).eq.'B') jqual = 5
      if (QUAL(1:1).eq.'C') jqual = 6
      if (QUAL(1:1).eq.'D') jqual = 7


      hqual = -1
      if (QUAL(2:2).eq.'X') hqual = 0
      if (QUAL(2:2).eq.'U') hqual = 1
      if (QUAL(2:2).eq.'F') hqual = 2
      if (QUAL(2:2).eq.'E') hqual = 3
      if (QUAL(2:2).eq.'A') hqual = 4
      if (QUAL(2:2).eq.'B') hqual = 5
      if (QUAL(2:2).eq.'C') hqual = 6
      if (QUAL(2:2).eq.'D') hqual = 7

      kqual = -1
      if (QUAL(3:3).eq.'X') kqual = 0
      if (QUAL(3:3).eq.'U') kqual = 1
      if (QUAL(3:3).eq.'F') kqual = 2
      if (QUAL(3:3).eq.'E') kqual = 3
      if (QUAL(3:3).eq.'A') kqual = 4
      if (QUAL(3:3).eq.'B') kqual = 5
      if (QUAL(3:3).eq.'C') kqual = 6
      if (QUAL(3:3).eq.'D') kqual = 7

      if (jqual.eq.-1.or.hqual.eq.-1.or.kqual.eq.-1) then
         write(*,*) 'jqual: ',jqual,QUAL(1:1)
         write(*,*) 'hqual: ',hqual,QUAL(2:2)
         write(*,*) 'kqual: ',kqual,QUAL(3:3)
         stop
         endif

      imj = abs(jmags)*8 + jqual
      imh = abs(hmags)*8 + hqual
      imk = abs(kmags)*8 + kqual

      if (jmags.lt.0) imj = -imj
      if (hmags.lt.0) imh = -imh
      if (kmags.lt.0) imk = -imk


c
c-------------------------------------------------------------
c
      do b = 1, 32
         b32_dei(b) = 0
         b32_rai(b) = 0
         enddo
      do b = 1, 16
         b16_mji(b) = 0
         b16_mhi(b) = 0
         b16_mki(b) = 0
         enddo

      call i4_to_b32(iras,b32_rai)
      call i4_to_b32(ides,b32_dei)


      call i2_to_b16(imj,b16_mji)
      call i2_to_b16(imh,b16_mhi)
      call i2_to_b16(imk,b16_mki)

c
c-------------------------------------------------------------
c
      call bbbbbbbb_to_b1(buff13(01),
     .        b32_dei(05),b32_dei(06),b32_dei(07),b32_dei(08),     ! 1st nybble in de_a 
     .        b32_dei(09),b32_dei(10),b32_dei(11),b32_dei(12))     ! 2nd nybble in de_a

      call bbbbbbbb_to_b1(buff13(02),
     .        b32_dei(13),b32_dei(14),b32_dei(15),b32_dei(16),     ! 3rd nybble in de_a
     .        b32_rai(05),b32_rai(06),b32_rai(07),b32_rai(08))     ! 1st nybble in ra_a  
 
      call bbbbbbbb_to_b1(buff13(03),
     .        b32_rai(09),b32_rai(10),b32_rai(11),b32_rai(12),     ! 2nd nybble in ra_a   
     .        b32_rai(13),b32_rai(14),b32_rai(15),b32_rai(16))     ! 3rd nybble in ra_a

      call bbbbbbbb_to_b1(buff13(04),                          
     .        b16_mji(01), b16_mhi(01), b16_mki(01),               ! 1st byte in mjhk
     .        b16_mji(02), b16_mhi(02), b16_mki(02),               ! 2nd byte in mjhk
     .        b16_mji(03), b16_mhi(03))                            ! 3rd byte in mjhk (2/3)

      call bbbbbbbb_to_b1(buff13(05),                          
     .                                  b16_mki(03),               ! 3rd byte in mjhk (1/3)
     .        b16_mji(04), b16_mhi(04), b16_mki(04),               ! 4th byte in mjhk
     .        b16_mji(05), b16_mhi(05), b16_mki(05),               ! 5th byte in mjhk
     .        b16_mji(06))                                         ! 6th byte in mjhk (1/3)

      call bbbbbbbb_to_b1(buff13(06),                          
     .                     b16_mhi(06), b16_mki(06),               ! 6th byte in mjhk (2/3)
     .        b16_mji(07), b16_mhi(07), b16_mki(07),               ! 7th byte in mjhk 
     .        b16_mji(08), b16_mhi(08), b16_mki(08))               ! 8th byte in mjhk

      call bbbbbbbb_to_b1(buff13(07),                          
     .        b16_mji(09), b16_mji(10), b16_mji(11), b16_mji(12),  ! last byte in mj
     .        b16_mji(13), b16_mji(14), b16_mji(15), b16_mji(16)) 

      call bbbbbbbb_to_b1(buff13(08),                          
     .        b16_mhi(09), b16_mhi(10), b16_mhi(11), b16_mhi(12),  ! last byte in mh
     .        b16_mhi(13), b16_mhi(14), b16_mhi(15), b16_mhi(16)) 
      
      call bbbbbbbb_to_b1(buff13(09),                          
     .        b16_mki(09), b16_mki(10), b16_mki(11), b16_mki(12),  ! last byte in mk
     .        b16_mki(13), b16_mki(14), b16_mki(15), b16_mki(16)) 

      call bbbbbbbb_to_b1(buff13(10),
     .        b32_rai(17), b32_rai(18), b32_rai(19), b32_rai(20), 
     .        b32_rai(21), b32_rai(22), b32_rai(23), b32_rai(24)) 

      call bbbbbbbb_to_b1(buff13(11),
     .        b32_rai(25), b32_rai(26), b32_rai(27), b32_rai(28), 
     .        b32_rai(29), b32_rai(30), b32_rai(31), b32_rai(32)) 

      call bbbbbbbb_to_b1(buff13(12),
     .        b32_dei(17), b32_dei(18), b32_dei(19), b32_dei(20), 
     .        b32_dei(21), b32_dei(22), b32_dei(23), b32_dei(24)) 

      call bbbbbbbb_to_b1(buff13(13),
     .        b32_dei(25), b32_dei(26), b32_dei(27), b32_dei(28), 
     .        b32_dei(29), b32_dei(30), b32_dei(31), b32_dei(32)) 

      return
      end



c---------------------------------------------------
c
c designed to compress better...
c
c
      subroutine buff2rdm(buff13,ira,ide,imj,imh,imk,aaa)
     .                   
      implicit none

      byte        buff13(13)
      integer*4   ira, ide
      integer*2   imj, imh, imk
      character*3 aaa

      integer*4 ira1, ira2, irau
      integer*4 ide1, ide2, ideu
      integer*2 imju, jqual
      integer*2 imhu, hqual
      integer*2 imku, kqual

      byte b

      integer*4 irax
      integer*4 idex

      integer*2 jmagu, jmagx
      integer*2 hmagu, hmagx
      integer*2 kmagu, kmagx

      byte   b32_rai(32)
      byte   b32_dei(32)
      byte   b16_mji(16)
      byte   b16_mhi(16)
      byte   b16_mki(16)

      real*8 ra, de
      real*8 jmag, hmag, kmag

c---------------------------------------------------------------------------
c
c now, back out what we put in, using only the buff information
c
      do b = 1, 32
         b32_dei(b) = 0
         b32_rai(b) = 0
         enddo
      do b = 1, 16
         b16_mji(b) = 0
         b16_mhi(b) = 0
         b16_mki(b) = 0
         enddo
      ira1 = 0
      ira2 = 0
      ide1 = 0
      ide2 = 0
      imj  = 0
      imh  = 0
      imk  = 0
      jmagu = 0
      hmagu = 0
      kmagu = 0
      AAA = '   '


c---------------------------------------------------------------------------
c
      call b1_to_bbbbbbbb(buff13(01),
     .        b32_dei(05),b32_dei(06),b32_dei(07),b32_dei(08),     ! 1st nybble in de_a 
     .        b32_dei(09),b32_dei(10),b32_dei(11),b32_dei(12))     ! 2nd nybble in de_a

      call b1_to_bbbbbbbb(buff13(02),
     .        b32_dei(13),b32_dei(14),b32_dei(15),b32_dei(16),     ! 3rd nybble in de_a
     .        b32_rai(05),b32_rai(06),b32_rai(07),b32_rai(08))     ! 1st nybble in ra_a  
 
      call b1_to_bbbbbbbb(buff13(03),
     .        b32_rai(09),b32_rai(10),b32_rai(11),b32_rai(12),     ! 2nd nybble in ra_a   
     .        b32_rai(13),b32_rai(14),b32_rai(15),b32_rai(16))     ! 3rd nybble in ra_a

      call b1_to_bbbbbbbb(buff13(04),                          
     .        b16_mji(01), b16_mhi(01), b16_mki(01),               ! 1st byte in mjhk
     .        b16_mji(02), b16_mhi(02), b16_mki(02),               ! 2nd byte in mjhk
     .        b16_mji(03), b16_mhi(03))                            ! 3rd byte in mjhk (2/3)

      call b1_to_bbbbbbbb(buff13(05),                          
     .                                  b16_mki(03),               ! 3rd byte in mjhk (1/3)
     .        b16_mji(04), b16_mhi(04), b16_mki(04),               ! 4th byte in mjhk
     .        b16_mji(05), b16_mhi(05), b16_mki(05),               ! 5th byte in mjhk
     .        b16_mji(06))                                         ! 6th byte in mjhk (1/3)

      call b1_to_bbbbbbbb(buff13(06),                          
     .                     b16_mhi(06), b16_mki(06),               ! 6th byte in mjhk (2/3)
     .        b16_mji(07), b16_mhi(07), b16_mki(07),               ! 7th byte in mjhk 
     .        b16_mji(08), b16_mhi(08), b16_mki(08))               ! 8th byte in mjhk

      call b1_to_bbbbbbbb(buff13(07),                          
     .        b16_mji(09), b16_mji(10), b16_mji(11), b16_mji(12),  ! last byte in mj
     .        b16_mji(13), b16_mji(14), b16_mji(15), b16_mji(16)) 

      call b1_to_bbbbbbbb(buff13(08),                          
     .        b16_mhi(09), b16_mhi(10), b16_mhi(11), b16_mhi(12),  ! last byte in mh
     .        b16_mhi(13), b16_mhi(14), b16_mhi(15), b16_mhi(16)) 
      
      call b1_to_bbbbbbbb(buff13(09),                          
     .        b16_mki(09), b16_mki(10), b16_mki(11), b16_mki(12),  ! last byte in mk
     .        b16_mki(13), b16_mki(14), b16_mki(15), b16_mki(16)) 

      call b1_to_bbbbbbbb(buff13(10),
     .        b32_rai(17), b32_rai(18), b32_rai(19), b32_rai(20), 
     .        b32_rai(21), b32_rai(22), b32_rai(23), b32_rai(24)) 

      call b1_to_bbbbbbbb(buff13(11),
     .        b32_rai(25), b32_rai(26), b32_rai(27), b32_rai(28), 
     .        b32_rai(29), b32_rai(30), b32_rai(31), b32_rai(32)) 

      call b1_to_bbbbbbbb(buff13(12),
     .        b32_dei(17), b32_dei(18), b32_dei(19), b32_dei(20), 
     .        b32_dei(21), b32_dei(22), b32_dei(23), b32_dei(24)) 

      call b1_to_bbbbbbbb(buff13(13),
     .        b32_dei(25), b32_dei(26), b32_dei(27), b32_dei(28), 
     .        b32_dei(29), b32_dei(30), b32_dei(31), b32_dei(32)) 

      call b32_to_i4(b32_rai,irau)
      call b32_to_i4(b32_dei,ideu)

      call b16_to_i2(b16_mji,imj)
      call b16_to_i2(b16_mhi,imh)
      call b16_to_i2(b16_mki,imk)

      ira1 = irau/65536
      ira2 = irau-irau/65536*65536

      ide1 = ideu/65536
      ide2 = ideu-ideu/65536*65536

      ira  = ira1*36000 + ira2
      ide  = ide1*36000 + ide2

      jqual = abs(imj)-abs(imj)/8*8
      hqual = abs(imh)-abs(imh)/8*8
      kqual = abs(imk)-abs(imk)/8*8

      jmagu = abs(imj)/8
      hmagu = abs(imh)/8
      kmagu = abs(imk)/8
      if (imj.lt.0) jmagu = -jmagu
      if (imh.lt.0) hmagu = -hmagu
      if (imk.lt.0) kmagu = -kmagu

      if (jqual.eq.0) AAA(1:1) = 'X'
      if (jqual.eq.1) AAA(1:1) = 'U'
      if (jqual.eq.2) AAA(1:1) = 'F'
      if (jqual.eq.3) AAA(1:1) = 'E'
      if (jqual.eq.4) AAA(1:1) = 'A'
      if (jqual.eq.5) AAA(1:1) = 'B'
      if (jqual.eq.6) AAA(1:1) = 'C'
      if (jqual.eq.7) AAA(1:1) = 'D'

      if (hqual.eq.0) AAA(2:2) = 'X'
      if (hqual.eq.1) AAA(2:2) = 'U'
      if (hqual.eq.2) AAA(2:2) = 'F'
      if (hqual.eq.3) AAA(2:2) = 'E'
      if (hqual.eq.4) AAA(2:2) = 'A'
      if (hqual.eq.5) AAA(2:2) = 'B'
      if (hqual.eq.6) AAA(2:2) = 'C'
      if (hqual.eq.7) AAA(2:2) = 'D'

      if (kqual.eq.0) AAA(3:3) = 'X'
      if (kqual.eq.1) AAA(3:3) = 'U'
      if (kqual.eq.2) AAA(3:3) = 'F'
      if (kqual.eq.3) AAA(3:3) = 'E'
      if (kqual.eq.4) AAA(3:3) = 'A'
      if (kqual.eq.5) AAA(3:3) = 'B'
      if (kqual.eq.6) AAA(3:3) = 'C'
      if (kqual.eq.7) AAA(3:3) = 'D'

      imj = jmagu
      imh = hmagu
      imk = kmagu

      return
      end


c----------------------------------------------
c
      subroutine i4_to_b32(i4,b32)
      implicit none
    
      integer*4 i4
      byte      b32(32)

      integer*4 i4u

      i4u = i4

      b32(01) = 0
      if (i4u.lt.0) then
         b32(01) = 0
         i4u = -i4u
         endif

      b32(02) = i4u/2**30
      b32(03) = (i4u-i4u/2**30*2**30)/2**29
      b32(04) = (i4u-i4u/2**29*2**29)/2**28
      b32(05) = (i4u-i4u/2**28*2**28)/2**27
      b32(06) = (i4u-i4u/2**27*2**27)/2**26
      b32(07) = (i4u-i4u/2**26*2**26)/2**25
      b32(08) = (i4u-i4u/2**25*2**25)/2**24
      b32(09) = (i4u-i4u/2**24*2**24)/2**23
      b32(10) = (i4u-i4u/2**23*2**23)/2**22
      b32(11) = (i4u-i4u/2**22*2**22)/2**21
      b32(12) = (i4u-i4u/2**21*2**21)/2**20
      b32(13) = (i4u-i4u/2**20*2**20)/2**19
      b32(14) = (i4u-i4u/2**19*2**19)/2**18
      b32(15) = (i4u-i4u/2**18*2**18)/2**17
      b32(16) = (i4u-i4u/2**17*2**17)/2**16
      b32(17) = (i4u-i4u/2**16*2**16)/2**15
      b32(18) = (i4u-i4u/2**15*2**15)/2**14
      b32(19) = (i4u-i4u/2**14*2**14)/2**13
      b32(20) = (i4u-i4u/2**13*2**13)/2**12
      b32(21) = (i4u-i4u/2**12*2**12)/2**11
      b32(22) = (i4u-i4u/2**11*2**11)/2**10
      b32(23) = (i4u-i4u/2**10*2**10)/2**09
      b32(24) = (i4u-i4u/2**09*2**09)/2**08
      b32(25) = (i4u-i4u/2**08*2**08)/2**07
      b32(26) = (i4u-i4u/2**07*2**07)/2**06
      b32(27) = (i4u-i4u/2**06*2**06)/2**05
      b32(28) = (i4u-i4u/2**05*2**05)/2**04
      b32(29) = (i4u-i4u/2**04*2**04)/2**03
      b32(30) = (i4u-i4u/2**03*2**03)/2**02
      b32(31) = (i4u-i4u/2**02*2**02)/2**01
      b32(32) = (i4u-i4u/2*2)

      return
      end




c----------------------------------------------
c
      subroutine i2_to_b16(i2,b16)
      implicit none
    
      integer*2 i2
      byte      b16(16)

      integer*2  i2u
      integer b


      i2u = i2
      b16(01) = 0
      if (i2u.lt.0) then
         b16(01) = 1
         i2u = -i2u
         endif

      b16(02) = (i2u-i2u/2**15*2**15)/2**14
      b16(03) = (i2u-i2u/2**14*2**14)/2**13
      b16(04) = (i2u-i2u/2**13*2**13)/2**12
      b16(05) = (i2u-i2u/2**12*2**12)/2**11
      b16(06) = (i2u-i2u/2**11*2**11)/2**10
      b16(07) = (i2u-i2u/2**10*2**10)/2**09
      b16(08) = (i2u-i2u/2**09*2**09)/2**08
      b16(09) = (i2u-i2u/2**08*2**08)/2**07
      b16(10) = (i2u-i2u/2**07*2**07)/2**06
      b16(11) = (i2u-i2u/2**06*2**06)/2**05
      b16(12) = (i2u-i2u/2**05*2**05)/2**04
      b16(13) = (i2u-i2u/2**04*2**04)/2**03
      b16(14) = (i2u-i2u/2**03*2**03)/2**02
      b16(15) = (i2u-i2u/2**02*2**02)/2**01
      b16(16) = (i2u-i2u/2*2)

      return
      end


c----------------------------------------------
c
      subroutine b1_to_b8(b1,b8)
      implicit none
    
      byte      b1
      byte      b8(8)

      integer   b1u

      b1u = b1

      if (b1u.lt.0) b1u = b1u + 256
 
      b8(01) = (b1u-b1u/2**08*2**08)/2**07
      b8(02) = (b1u-b1u/2**07*2**07)/2**06
      b8(03) = (b1u-b1u/2**06*2**06)/2**05
      b8(04) = (b1u-b1u/2**05*2**05)/2**04
      b8(05) = (b1u-b1u/2**04*2**04)/2**03
      b8(06) = (b1u-b1u/2**03*2**03)/2**02
      b8(07) = (b1u-b1u/2**02*2**02)/2**01
      b8(08) = (b1u-b1u/2**01*2**01)/2**00

      return
      end


c----------------------------------------------
c
      subroutine b8_to_b1(b8,b1)
      implicit none
    
      byte      b8(8)
      byte      b1


      b1 = b8(8)*2**0
     .   + b8(7)*2**1
     .   + b8(6)*2**2
     .   + b8(5)*2**3
     .   + b8(4)*2**4
     .   + b8(3)*2**5    
     .   + b8(2)*2**6

      if (b8(1).eq.1) b1 = -b1

      return
      end


c----------------------------------------------
c
      subroutine b16_to_i2(b16,i2)
      implicit none
    
      byte      b16(16)
      integer*2 i2


      i2 = b16(16)*2**00
     .   + b16(15)*2**01
     .   + b16(14)*2**02
     .   + b16(13)*2**03
     .   + b16(12)*2**04
     .   + b16(11)*2**05    
     .   + b16(10)*2**06
     .   + b16(09)*2**07     
     .   + b16(08)*2**08     
     .   + b16(07)*2**09
     .   + b16(06)*2**10     
     .   + b16(05)*2**11     
     .   + b16(04)*2**12     
     .   + b16(03)*2**13     
     .   + b16(02)*2**14     

      if (b16(01).eq.1) i2 = -i2
 
      return
      end


c----------------------------------------------
c
      subroutine b32_to_i4(b32,i4)
      implicit none
    
      byte      b32(32)
      integer*4 i4

      i4 = b32(32)*2**00
     .   + b32(31)*2**01
     .   + b32(30)*2**02
     .   + b32(29)*2**03
     .   + b32(28)*2**04
     .   + b32(27)*2**05    
     .   + b32(26)*2**06
     .   + b32(25)*2**07     
     .   + b32(24)*2**08     
     .   + b32(23)*2**09
     .   + b32(22)*2**10     
     .   + b32(21)*2**11     
     .   + b32(20)*2**12     
     .   + b32(19)*2**13     
     .   + b32(18)*2**14     
     .   + b32(17)*2**15     
     .   + b32(16)*2**16    
     .   + b32(15)*2**17    
     .   + b32(14)*2**18    
     .   + b32(13)*2**19     
     .   + b32(12)*2**20     
     .   + b32(11)*2**21   
     .   + b32(10)*2**22    
     .   + b32(09)*2**23     
     .   + b32(08)*2**24     
     .   + b32(07)*2**25     
     .   + b32(06)*2**26     
     .   + b32(05)*2**27     
     .   + b32(04)*2**28     
     .   + b32(03)*2**29     
     .   + b32(02)*2**30     
      if (b32(01).eq.1) i4 = -i4

      return
      end


c----------------------------------------------------------
c
      subroutine bbbbbbbb_to_b1(b1,ba,bb,bc,bd,be,bf,bg,bh)
      implicit none

      byte b1
      byte ba,bb,bc,bd
      byte be,bf,bg,bh

      integer b1u
     
      byte bax,bbx,bcx,bdx
      byte bex,bfx,bgx,bhx

      b1u = 0
      b1u = b1u 
     .    + bh 
     .    + bg*002  
     .    + bf*004
     .    + be*008
     .    + bd*016
     .    + bc*032  
     .    + bb*064  
     .    + ba*128

      b1 = b1u 

      call b1_to_bbbbbbbb(b1,bax,bbx,bcx,bdx,bex,bfx,bgx,bhx)

      if (ba.ne.bax.or.
     .    bb.ne.bbx.or.
     .    bc.ne.bcx.or.
     .    bd.ne.bdx.or.
     .    be.ne.bex.or.
     .    bf.ne.bfx.or.
     .    bg.ne.bgx.or.
     .    bh.ne.bhx) then
             print*,'b1 : ',b1
             print*,'b1u: ',b1u
             print*,'  a: ',ba,bax
             print*,'  b: ',bb,bbx
             print*,'  c: ',bc,bcx
             print*,'  d: ',bd,bdx
             print*,'  e: ',be,bex
             print*,'  f: ',bf,bfx
             print*,'  g: ',bg,bgx
             print*,'  h: ',bh,bhx
             stop
             endif

      return
      end

c----------------------------------------------------------
c
      subroutine b1_to_bbbbbbbb(b1,ba,bb,bc,bd,be,bf,bg,bh)
      implicit none

      byte b1
      byte ba,bb,bc,bd
      byte be,bf,bg,bh

      integer*2 b1u

      b1u = b1
      if (b1u.lt.0) b1u = b1u + 256
 
      ba = (b1u-b1u/2**08*2**08)/2**07
      bb = (b1u-b1u/2**07*2**07)/2**06
      bc = (b1u-b1u/2**06*2**06)/2**05
      bd = (b1u-b1u/2**05*2**05)/2**04
      be = (b1u-b1u/2**04*2**04)/2**03
      bf = (b1u-b1u/2**03*2**03)/2**02
      bg = (b1u-b1u/2**02*2**02)/2**01
      bh = (b1u-b1u/2**01*2**01)/2**00
 
      return
      end


c-------------------------------------------
c
c this routine will convert a string into RA,
c it can handle both decimal and sexigesimal
c
      real*8 function string2ra(string)
      implicit none
      character*15 string

      logical iscolon
      character*15 stringc
      real rah, ram, ras 
      integer i

      iscolon = .false. 
      do i = 01, 15
         stringc(i:i) = string(i:i)
         if (stringc(i:i).eq.':') then
            iscolon = .true. 
            stringc(i:i) = ' '
            endif
         enddo

      if (.not.iscolon) then
         read(stringc,*) string2ra
         return
         endif
      if (iscolon) then
         read(stringc,*) RAH, RAM, RAS
         string2ra = (rah + ram/60. + ras/60./60.)*15.0
         return
         endif
          
      return
      end


c-------------------------------------------
c
c this routine will convert a string into RA,
c it can handle both decimal and sexigesimal
c
      real*8 function string2de(string)
      implicit none
      character*15 string

      logical iscolon
      logical isneg
      character*15 stringc
      real DED, DEM, DES
      integer i


      iscolon = .false. 
      isneg   = .false.
      do i = 01, 15
         stringc(i:i) = string(i:i)
         if (stringc(i:i).eq.':') then
            iscolon = .true. 
            stringc(i:i) = ' '
            endif
         if (stringc(i:i).eq.'-') isneg = .true.
         enddo

      if (.not.iscolon) then
         read(stringc,*) string2de
         return
         endif

      if (iscolon) then
         read(stringc,*) DED, DEM, DES
         string2de = (abs(DED) + DEM/60. + DES/60./60.)
         if (isneg) string2de = -string2de
         return
         endif
          
      return
      end



c-------------------------------------------
c
c this routine will convert a string into RA,
c it can handle both decimal and sexigesimal
c
      character*11 function ra2sexig(ra)
      implicit none
      real*8 ra

      integer rah, ram, rasi, rasf_x10
      real    ras, rasf

      rah = int(ra/15)
      ram = int((ra/15-rah)*60)
      ras = int((ra/15-rah-ram/60.)*60*60)
      rasi = ras
      rasf = ras-rasi
      rasf_x10 = rasf*100

      write(ra2sexig,'(i2.2,'':'',i2.2,'':'',i2.2,''.'',i2.2)')
     .    rah, ram, rasi, rasf_x10

      return
      end


c-------------------------------------------
c
c this routine will convert a string into RA,
c it can handle both decimal and sexigesimal
c
      character*12 function de2sexig(de)
      implicit none
      real*8 de

      real    deu
      integer ded, dem, desi, desf_x10
      real    des, desf
      character neg
 
      if (de.gt.0) neg = '+'
      if (de.lt.0) neg = '-'
      deu = abs(de)
      ded = int(deu)
      dem = int((deu-ded)*60)
      des = int((deu-ded-dem/60.)*60*60)
      desi = des
      desf = des-desi
      desf_x10 = desf*100

      write(de2sexig,'(a1,i2.2,'':'',i2.2,'':'',i2.2,''.'',i2.2)')
     .    neg, ded, dem, desi, desf_x10

      return
      end





      subroutine writfits_b1(FILE,pix,PXDIMX,PXDIMY)
      implicit none

      character*(*) FILE
      integer PXDIMX,PXDIMY
      byte pix(PXDIMX,PXDIMY)

      integer nbyte0
      integer nbyteE
      integer nbyte1
      integer nbyte2
      integer nbper
      integer i,ios

      character*2880 buffc
      byte buffb(2880)
      equivalence (buffc,buffb)

      integer ifirst, i1, i2

      integer np1, np2, npt

      character*70 HDR(25)
      common/HDR/HDR

      character*80 FILEU

      FILEU = FILE
      do i = 75,2,-1
         if (FILE(i:i+4).eq.'.fits') FILEU = FILE(1:i+4)
         enddo



c     write(*,'(''ENTER   writfits_b1: '',80a)') FILEU

      open(10,file=FILEU,
     .     status='unknown',
     .     err   =900,
     .     recl  =2880,
     .     form  ='UNFORMATTED',
     .     access='DIRECT')


      i = 1
      write(buffc( 0*80+1: 1*80),'(''SIMPLE  = T    '')')
      write(buffc( 1*80+1: 2*80),'(''BITPIX  =  8  '')')
      write(buffc( 2*80+1: 3*80),'(''NAXIS   = '',i12)') 2
      write(buffc( 3*80+1: 4*80),'(''NAXIS1  = '',i12)') PXDIMX
      write(buffc( 4*80+1: 5*80),'(''NAXIS2  = '',i12)') PXDIMY
      write(buffc( 5*80+1: 6*80),'(''DATATYPE= '',9a)') 
     .                          ' ''BYTE    '''
      write(buffc( 6*80+1: 7*80),'(''DATE    = '',11a)') 
     .                          ' ''28/01/00'''
      write(buffc( 7*80+1: 8*80),'(''BSCALE  = '',i12)') 00001
      write(buffc( 8*80+1: 9*80),'(''BZERO   = '',i12)') 00000
      write(buffc(09*80+1:10*80),'(''CRPIX1  = '',a20)') HDR(01)
      write(buffc(10*80+1:11*80),'(''CRPIX2  = '',a20)') HDR(02)
      write(buffc(11*80+1:12*80),'(''CRVAL1  = '',a20)') HDR(03)
      write(buffc(12*80+1:13*80),'(''CRVAL2  = '',a20)') HDR(04)
      write(buffc(13*80+1:14*80),'(''CTYPE1  = '',a20)') HDR(05)
      write(buffc(14*80+1:15*80),'(''CTYPE2  = '',a20)') HDR(06)
      write(buffc(15*80+1:16*80),'(''CD1_1   = '',a20)') HDR(07)
      write(buffc(16*80+1:17*80),'(''CD1_2   = '',a20)') HDR(08)
      write(buffc(17*80+1:18*80),'(''CD2_1   = '',a20)') HDR(09)
      write(buffc(18*80+1:19*80),'(''CD2_2   = '',a20)') HDR(10)
      write(buffc(19*80+1:20*80),'(''ORIENTAT= '',a20)') HDR(11)
      write(buffc(20*80+1:21*80),'(''PA_APER = '',a20)') HDR(12)
      write(buffc(21*80+1:22*80),'(''PA_V3   = '',a20)') HDR(13)
      write(buffc(22*80+1:23*80),'(''DATE-OBS= '',a20)') HDR(14)
      write(buffc(23*80+1:24*80),'(''TIME-OBS= '',a20)') HDR(15)
      write(buffc(24*80+1:25*80),'(''EXPTIME = '',a20)') HDR(16)
      write(buffc(25*80+1:26*80),'(''ROOTNAME= '',a20)') HDR(17)
      write(buffc(26*80+1:27*80),'(''TARGNAME= '',a20)') HDR(18)
      write(buffc(27*80+1:28*80),'(''RA_TARG = '',a20)') HDR(19)
      write(buffc(28*80+1:29*80),'(''DEC_TARG= '',a20)') HDR(20)
      write(buffc(29*80+1:30*80),'(''PROPOSID= '',a20)') HDR(21)
      write(buffc(30*80+1:31*80),'(''FILTER1 = '',a20)') HDR(22)
      write(buffc(31*80+1:32*80),'(''FILTER2 = '',a20)') HDR(23)
      write(buffc(32*80+1:33*80),'(''VAFACTOR= '',a20)') HDR(24)
      write(buffc(33*80+1:34*80),'(''COMMENT   '')')
      write(buffc(34*80+1:35*80),'(''COMMENT   '')')
      write(buffc(35*80+1:36*80),'(''END       '')')

c     print*,buffc
      write(10,rec=i,iostat=ios) buffc

      ifirst = i+1
      i1 = i
      i2 = i

      nbper  =   PXDIMX*PXDIMY
      npt    =   PXDIMX*PXDIMY
      nbyte1 = 1 
      nbyte2 = nbper
      i1 = i+1 + nbyte1/2880
      i2 = i+1 + nbyte2/2880


c      print*,'nbyte1: ',nbyte1
c      print*,'nbyte2: ',nbyte1
c      print*,'    i1: ',i1 
c      print*,'    i2: ',i2
      do i = i1, i2, 1
         nbyte0 = (i-ifirst)*2880+   1
         nbyteE = (i-ifirst)*2880+2880
         np1 = (nbyte0-nbyte1) + 1
         np2 = (nbyteE-nbyte1)+ 1
        !if (i.lt.0010) print*,'i: ',i,i1,i2,np1,np2,nbyte0
        !if (i.gt.2900) print*,'i: ',i,i1,i2,np1,np2,nbyte0
         call pix2buff_b1(pix,buffb,np1,npt)
         write(10,rec=i,iostat=ios) buffc
         enddo 

      close(10)

      return

 900  continue
      print*,'WRITFITS.f ERROR'
      stop

      end

      subroutine pix2buff_b1(pix,buffb,n1,nt)
      implicit none 
      byte pix(1)
      byte buffb(2880)

      integer   n1,nt
 
      integer i
      integer npu

      integer ipval

      do i = 1, 2880
         npu = (1+(n1-1)*2)+i-1
         npu = (1+(n1-1)*1)+i-1
         if (npu.ge.1.and.npu.le.nt) then
            ipval = pix(npu)
            buffb(i) = ipval
            endif
         enddo
 
      return
      end


c#include "/home/jayander/FORTRAN/MAC_OSX/BUFF/buff_work.f"
c#include "/home/jayander/FORTRAN/MAC_OSX/FITS/writfits_i4.f"
c#include "/home/jayander/FORTRAN/MAC_OSX/FITS/writfits_b1.f"
c#include "/home/jayander/FORTRAN/MAC_OSX/FITS/readfits_i4.f"


