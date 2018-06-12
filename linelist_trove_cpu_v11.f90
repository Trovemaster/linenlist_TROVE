program linelist_trove
 implicit none

 integer isym,j,Jf_,Ji_,ilevel,ilevelf_,ileveli_,igamma,igammaf_,igammai_,info,j0,i,k,jgamma
 integer nfiles, jmax, ilines
 integer nlevels,v_(30),itermi_,itermf_,igammav,sym_t,J_t,imis,v_norm(30),igammar,ipolyad,npolyad(0:100)

 real(8),allocatable :: energies(:),coeff_(:)
 integer,allocatable :: sym(:),Jktau(:,:),NN(:,:,:),Nlines(:),npolyadJG(:,:,:),numj(:,:)
 integer(8)  :: nswap,iswap,nswap_,iselect,nselect,old_swaps,kswaps
 real(8),allocatable :: acoef(:),tranfreq(:)
 integer(4),allocatable :: ilevelf(:),ileveli(:),igammaf(:),igammai(:),jf(:),ji(:),itermi(:),itermf(:)
 integer(4),allocatable :: ilevel_(:),v_l(:,:),v_norm_(:,:),igammav_(:),igammar_(:)
 character(4),allocatable :: gamma_(:),gammar_(:),gammav_(:)

 logical :: energeyfile_do = .false. , ifopen, ifappend = .false. , ifexist

 real(8) :: acoef_,abscoef,energy,energyf,energyi,tranfreq_,linestr,ZPE,coeff,nu1,nu2

 character(4) gamma,gammaf,gammai,gammav,jchar,gchar,gammar
 character(100) intfilename(4000),enrfilename,enrfilenameout,intfilenameout,dscrfilename,filename
 character(500) fname
 character(1)  :: ch_t1
 character(2)  :: ch_t2
 character(4)  :: fchar
 character(7)  :: fchar1,fchar2
 character(200) :: ch_t200
 character(400) :: form
 integer        :: unit_f,ilevelmax,nmodes,nsym
 character(len=4),allocatable   :: symm(:)
 integer ,allocatable     :: gns(:)

 type numT
      !
      integer,pointer   :: nn(:)
      integer,pointer   :: sym(:)
      integer,pointer   :: J(:)
      real(8),pointer   :: energy(:)
      !
 end type numT

 type(numT),allocatable :: level(:,:)

 real(8),parameter :: planck=6.6260693d-27,avogno=6.0221415d+23,vellgt=2.99792458d+10,boltz=1.380658d-16

    !read number of intensity files
    read*,nfiles
    !
    if (nfiles>4000) then 
      print('("Too many files (>400):",i)'),nfiles
      stop 'Too many files'
    endif
    !
    !read intensities filenames
    do i = 1,nfiles
     read*,intfilename(i)
    enddo
    !
    !read energies filename
    read*,enrfilename
    !
    !read descr-filename - in
     read*,dscrfilename
    !
    !read intensity filename - out
     read*,intfilenameout
    !
    !read energies filename - out 
     read*,enrfilenameout

    !read nmodes, Nsym
     read*,nmodes,Nsym

    allocate(symm(nsym),gns(nsym),stat=info); 
    if (info/=0) stop 'error: symm,gns  out of memory'

    !read nuclear spin statistics weights
    read*,(symm(isym),gns(isym),isym=1,nsym)

    !read the frequency window
    read*,nu1,nu2
    !
    print*,'range =  ',nu1,' - ',nu2

    !read the number of lines to keep in memory
    read*,nswap
    !
    print*,'nswap = ',nswap
    !
    !read*,old_swaps
    !
    !print*,"write after ",old_swaps,"swaps"
    !  
    !compute constants
    !ln2=log(2.0)
    !pi=acos(-1.0)
    !beta=planck*vellgt/(boltz*temp)
    !cmcoef=avogno/(8.0*pi*vellgt)
    !dpwcoef=sqrt(2.0*ln2*boltz*avogno)/vellgt
    !
    !read the energy file 
    open(unit=11,file=trim(enrfilename))
    !
    allocate(numj(0:200,nsym),stat=info); 
    if (info/=0) stop 'error: numj is out of memory'
    !
    i = 0
    numj = 0
    jmax = 0
    !
    print("('Counting states...')")
    !
    do
       !
       read(11,"(2x,a4,i7,f14.6,9x,i3)",end=11) gamma,ilevel,energy,j
! A1        1      0.000000   ( A1 ;  0  0  0) ( A1 ;   0   0   0)      1.00 (   0   0   0   0) (     1)

       !
       i = i + 1
       !
       do isym = 1,nsym
         if (trim( adjustl(gamma)) == symm(isym)) then 
           igamma = isym
           exit
         endif
       enddo
       !
       numj(j,igamma) = numj(j,igamma) + 1
       !
       jmax = max(jmax,j)
       !
       cycle
11     exit
    end do
    !
    print("(' jmax = ',i4)"),jmax
    !
    nlevels = i
    print("(' nlevels = ',i)"),nlevels
    !
    rewind(11)
    !
    allocate(sym(nlevels),energies(nlevels),Jktau(nlevels,3),NN(0:Jmax,nsym,2),npolyadJG(0:Jmax,nsym,0:100),stat=info); 
    if (info/=0) stop 'error: sym,v,n,l,Jktau is out of memory'
    !
    if ( all(trim(enrfilenameout)/=(/'none','NONE','null','NULL'/)) ) then  
      !
      energeyfile_do = .true.
      !
      open(unit=12,file=trim(enrfilenameout),action='write',status='replace')
      !
    endif
    !
    sym = 0
    NN(:,:,1) = 0
    NN(:,:,2) = nlevels
    !
    NN(0,1,1) = 1
    !
    sym_t = 1
    J_t =  0
    !
    print("('Reading energies...')")
    !
    ilevelmax = 0
    npolyad = 0
    npolyadJG = 0

    allocate(gamma_(nlevels),ilevel_(nlevels),gammar_(nlevels),gammav_(nlevels),v_l(nmodes,nlevels),coeff_(nlevels),v_norm_(nmodes,nlevels),igammav_(nlevels),igammar_(nlevels),stat=info)
    if (info/=0) stop 'error: swap-objects  out of memory'
    !
    if (jmax>200) stop 'Jmax>200'
       !print*,gamma,ilevel,energies(i),gammar,Jktau(i,1:3)
       !
!  A1        1      0.000000   ( A1 ;  0  0  0) ( A1 ;   0   0   0   0   0   0)      0.99 (   0   0   0   0   0   0   0) (     1)
!  A1        1      0.000000   ( A1 ;  0  0  0 ) ( A1 ;   0   0   0   0   0   0   0   0   0 )      0.98 (   0   0   0   0   0   0   0   0   0   0 ) (     1 )
    !
    do i = 1,nlevels
       !
       read(11,"(a200)") ch_t200
       !
       read(ch_t200,"(2x,a4,i7,f14.6,4x,a4,1x,3i3,4x,a4,1x,<nmodes>i4,7x,f5.2,2x,<nmodes>i4)") gamma,ilevel,energies(i),gammar,Jktau(i,1:3),gammav,v_(1:nmodes),coeff,v_norm(1:nmodes)
       !
       gamma_(i)       = gamma
       ilevel_(i)      = ilevel
       gammar_(i)      = gammar
       gammav_(i)      = gammav
       v_l(1:nmodes,i) = v_(1:nmodes)
       coeff_(i)       = coeff
       v_norm_(1:nmodes,i)      = v_norm(1:nmodes)
       !
    enddo
    !
    print("('Converting symmetries...')")
    !
    do i = 1,nlevels
       !
       gamma        = gamma_(i) 
       ilevel       = ilevel_(i)
       gammar       = gammar_(i)
       gammav       = gammav_(i)
       v_(1:nmodes) = v_l(1:nmodes,i)
       coeff        = coeff_(i) 
       v_norm(1:nmodes) = v_norm_(1:nmodes,i)
       !
       v_norm(1:nmodes) = v_(1:nmodes)
       !
       !print("(3i6)"),v_norm
       !
       !
       do j = 1,nsym
         if (trim( adjustl(gamma)) == symm(j)) then 
           sym(i) = j
           exit
         endif
       enddo
       !
       if (sym(i)==0) then
          write(6,"('wrong total symmetry of the level: ',a,i6)") adjustl(gamma),sym(i)
          stop 'wrong total symmetry of a level'
       endif
       !
       ! polyad number
       !
       igamma = sym(i)
       Ji_     = Jktau(i,1)
       !
       ipolyad = 2*(v_norm(1)+v_norm(3))+v_norm(2)
       !
       npolyad(ipolyad) = npolyad(ipolyad) + 1
       npolyadJG(Ji_,igamma,ipolyad) = npolyadJG(Ji_,igamma,ipolyad) + 1
       !
       if (sym(i)==0) stop 'wrong symmetry of the level'
       !
       if (sym_t /= igamma.or.J_t/=Ji_) then 
          NN(Ji_,igamma,1) = i
          NN(J_t,sym_t,2) = i-1
          !
          sym_t = igamma
          J_t = Ji_
          !
       endif
       !
       !if (jktau(i,1)>39) cycle
       !
       ilevelmax = max(ilevel,ilevelmax)
       !
       igammav = 0
       !
       do j = 1,nsym
         if (trim( adjustl(gammav)) == symm(j)) then 
           igammav = j
           exit
         endif
       enddo
       !
       igammav_(i) = igammav
       !
       if (igammav==0) then
          write(6,"('wrong vib symmetry of the level: ',a,i6)") adjustl(gammav),j
          stop 'wrong vib symmetry of the level'
       endif
       !
       igammar = 0
       !
       do j = 1,nsym
         if (trim( adjustl(gammar)) == symm(j)) then 
           igammar = j
           exit
         endif
       enddo
       !
       igammar_(i) = igammar
       !
       if (igammar==0) then 
           print*,'wrong rot symmetry of the level',gammar
           write(6,"(i12,f12.6,i6,i7,i4,2x,<nmodes>i4,3x,4i4,3x,<nmodes>i4,2x,i4)") & 
                 i,energies(i),gns(sym(i))*(2*jktau(i,1)+1),jktau(i,1),sym(i),v_norm(1:nmodes),igammar,jktau(i,1:3),v_(1:nmodes),igammav
         stop 'wrong rot symmetry of the level'
       endif
       ! 
    enddo
    !
    print("('Printing states file...')")
    !
    if (energeyfile_do) then
       !
       do i = 1,nlevels
         !
         !if (energies(i)<=99999.99) then 
            !
            write(12,"(i12,1x,f12.6,1x,i6,1x,i7,1x,i4,2x,<nmodes>i4,1x,i4,3x,4i4,2x,i8,2x,f5.2,3x,<nmodes>i4)") & 
                       !i,energies(i),gns(sym(i))*(2*jktau(i,1)+1),jktau(i,1),sym(i),v_norm(1:4),igammav,jktau(i,1:3),igammar,ilevel,coeff,v_(1:4)
                       !
                       i,energies(i),gns(sym(i))*(2*jktau(i,1)+1),jktau(i,1),sym(i),v_norm_(1:nmodes,i),igammav_(i),jktau(i,1:3),igammar_(i),0,coeff_(i),v_l(1:nmodes,i)
       
         !else
         !   write(12,"(i12,1x,e12.6,1x,i6,1x,i7,1x,i4,2x,<nmodes>i4,1x,i4,3x,4i4,2x,i8,2x,f5.2,3x,<nmodes>i4)") & 
         !              !i,energies(i),gns(sym(i))*(2*jktau(i,1)+1),jktau(i,1),sym(i),v_norm(1:4),igammav,jktau(i,1:3),igammar,ilevel,coeff,v_(1:4)
         !              !
         !              i,energies(i),gns(sym(i))*(2*jktau(i,1)+1),jktau(i,1),sym(i),v_norm_(1:nmodes,i),igammav_(i),jktau(i,1:3),igammar_(i),0,coeff_(i),v_l(1:nmodes,i)
         !       
         !endif
         ! 
       end do
       !
    endif
    !
    deallocate(gamma_,ilevel_,gammar_,gammav_,v_l,coeff_,v_norm_,igammav_,igammar_,stat=info)
    if (info/=0) stop 'error: swap-objects  out of memory'

    !
    print("('Maximal number of levels in a block = ',i9)"), ilevelmax
    !
    NN(J_t,sym_t,2) = nlevels
    !
    print*,"jmax= ",jmax
    !
    print("(' Energy ranges:')")
    !
    do j = 0,jmax
      !
      do igamma = 1,nsym
        if ( NN(j,igamma,1)/=0 ) then
          if (igamma/=nsym) then
             print("(i4,1x,i3,2x,'[',i8,'...',i8,'] = ',i8)"),j,igamma,NN(j,igamma,1:2),NN(j,igamma,2)-NN(j,igamma,1)+1
          else
             print("(i4,1x,i3,2x,'[',i8,'...',i8,'] = ',i8,' tot = ',i8)"),j,igamma,NN(j,igamma,1:2),NN(j,igamma,2)-NN(j,igamma,1)+1,NN(j,igamma,2)-NN(j,1,1)+1
          endif
          !
        endif
      enddo
      !
    enddo
    !
    close(11)
    !
    if (energeyfile_do) close(12)
    !
    !deallocate(nvib,lvib,tauvib,vvib,symvib)
    !
    print("(' N of levels =  ',i)"),nlevels
    !
    print("('Read in the descr-files and re-assign the numbering of the levels ..')")
    !
    ! finish if the intensity output is not requested
    !
    if (trim(intfilenameout)=="none".or.trim(intfilenameout)=="NONE") stop 
    !
    allocate(level(0:jmax,nsym))
    !
    do j = 0,jmax
      !
      do jgamma = 1,nsym
        !
        if (abs(gns(jgamma))<=1e-3.and.j/=0) cycle
        !
        allocate(level(j,jgamma)%nn(numj(j,jgamma)),level(j,jgamma)%energy(numj(j,jgamma)),stat=info); 
        if (info/=0) stop 'error: level%nn is out of memory'
        !
        level(j,jgamma)%nn = 0
        level(j,jgamma)%energy = 0
        !
      enddo
      !
    enddo
    k = 0
    ZPE = 0
    !
    do j = 0,jmax
      !
      write(jchar, '(i4)') j
      !
      do jgamma = 1,nsym
        !
        if (abs(gns(jgamma))<=1e-3.and.j/=0) cycle
        !
        write(gchar, '(i2)') jgamma
        !
        filename = trim(dscrfilename)//trim(adjustl(jchar))//"_"//trim(adjustl(gchar))//'.chk'
        !
        print("( 'process ',a)"), filename
        !
        open(unit=14,file=trim(filename),action='read',status='old')
        !
        do j0 = 1,14+nmodes
          read(14,"(a400)") form
        enddo
        !
        i = 0
        !
        do 
           !
           read(14,"(a400)",end=14) form
           !
           if (form(1:3)=='End') exit
           !
           read(form,*) ilevel,igamma,ilevel,igammai_,energy
           !
           !print*,ilevel,igamma,ilevel,igammai_,energy
           !
           if (igamma/=jgamma) then
             print("( 'illegal symmetry ',i8,' from the file with symmetry ',i,' = ',a  )"), igamma,jgamma,trim(filename)
             stop 'illegal symmetry'
           endif 
           !
           if (j==0.and.jgamma==1.and.ilevel==1) then
             ZPE = energy
             print*,"ZPE=",ZPE
           endif
           !
           i = i + 1
           !
           do while(k<nlevels)
             !
             k  = k + 1 
             !
             !print*,jktau(k,1),sym(k),energies(k)
             !
             if (j==jktau(k,1).and.igamma == sym(k) .and. abs(energy-ZPE-energies(k))<=2e-4) then
               !
               level(j,igamma)%nn(ilevel) = k
               level(j,igamma)%energy(ilevel) = energy-ZPE
               !
               exit
               !
             endif
             !
           enddo
           !
           if (level(j,igamma)%nn(ilevel)==0) then 
             !
             print("( 'cannot find a match from the descr-file for ilevel=',i8,' energy =',f14.4 )"), ilevel,energy-ZPE
             write(6,"('j,igamma = ',2i9)") j,gamma
             stop 'cannot find a match for the descr-file'
             !
           endif
           !
           cycle
        14 exit
        end do
        !
        close(14)
        !
      enddo
      !
    enddo
    !
    !deallocate(sym,energies,Jktau)
    !
    print("('Generate the Transition file ...')")
    !
    !open(unit=13,file=trim(intfilenameout),action='write',status='replace')
    !
    !start loop over all transitions
    !
    allocate(Nlines(0:jmax))
    !
    allocate(tranfreq(nswap),ilevelf(nswap),jf(nswap),igammaf(nswap),ileveli(nswap),ji(nswap),igammai(nswap),acoef(nswap),itermi(nswap),itermf(nswap),stat=info)
    if (info/=0) stop 'error: swap-objects  out of memory'
    !
    ilines = 0
    Nlines = 0
    !
    write(fchar1, '(f7.0)') nu1
    write(fchar2, '(f7.0)') nu2
    !write(fchar, '(i4)') old_swaps
    !
    fname = trim(intfilenameout)//'_'//trim(adjustl(fchar1))//'-'//trim(adjustl(fchar2))//'_'//trim(adjustl(fchar))//'.trans'
    !
    INQUIRE (FILE = fname, opened = ifopen)
    INQUIRE (FILE = fname, exist = ifexist)
    !
    unit_f = 20
    !
    !if (.not.ifopen)  open(unit=unit_f,file=fname,action='write',buffered='yes',status='replace')
    !
    !open(unit=unit_f,file=fname,action='write',buffered='yes',status='replace')
    !
    do i = 1,nfiles
      !
      open(unit=1,file=trim(intfilename(i)))
      !
      print("( 'process ',a,3x,i)"), intfilename(i),ilines
      !
      fname = trim(intfilename(i))//'_'//trim(adjustl(fchar1))//'-'//trim(adjustl(fchar2))//'.trans'
      !
      unit_f = 20
      !
      open(unit=unit_f,file=fname,action='write',buffered='yes',status='replace')
      !
      ilines = 0
      imis = 0
      kswaps = 0 
      !
      loop_file : do
        !
        nswap_ = nswap
        iselect = 0
        nselect = 0
        !
        kswaps = kswaps + 1
        !
        print*,"Reading transitions ...  ",kswaps
        !
        do iswap = 1,nswap
          !
          !read(1,*,end=20) ilevelf(iswap),ileveli(iswap),acoef(iswap)
     
          !   read new line from intensities file
          !
          iselect = iselect + 1
          !
          read(1,*,end=20) tranfreq(iselect),ilevelf(iselect),jf(iselect),igammaf(iselect),ch_t2,ileveli(iselect),ji(iselect),igammai(iselect),acoef(iselect)
          !
          j = max(Ji(iselect),Jf(iselect))
          !
          !write(6,"(i8,f16.6,7i7)") iselect,tranfreq(iselect),ilevelf(iselect),jf(iselect),igammaf(iselect),ileveli(iselect),ji(iselect),igammai(iselect)
          !
          if ( tranfreq(iselect)<nu1-1e-6 .or. tranfreq(iselect)>=(nu2-1e-6) .or. j>jmax) iselect = iselect - 1
          !
          !if (j>jmax) cycle
          !
          Nlines(j)= Nlines(j) + 1
          !
          !read(1,*,end=20) & 
          !         jf,igammaf,ji,igammai,energyf,energyi,tranfreq,acoef,abscoef,ilevelf,ileveli
     
          !      1  2            0  1       2581.0051       -0.0000   2581.0051     7.26055990E-02   3.13605334E+02       19        1 ||
          ! 
          !read(1,*,end=20) jf,ilevelf,gammaf,ji,ileveli,gammai,tranfreq,acoef,linestr,abscoef,energyf,energyi
     
          !
          cycle 
          !
        20 continue 
          !
          nswap_ = iswap-1
          iselect = iselect - 1
          !
          exit
          !
        enddo
        !
        nselect = iselect
        !
        if (nswap_<1) then 
           write(6,"('nswap = ')") nswap_
           stop 'nswap_=0'
        endif
        !
        !!!!if ( kswaps<=old_swaps ) cycle
        !
        print*,"Correlate states IDs ..."
        !
        if (nswap_<nswap) print('(a,i14,a)'),"last swap with ... ",nswap_," lines"
        !
        !print('(a,i14)'),"Correlate states IDs for ",nselect," transitions ..."
        !
        !$omp parallel do private(iswap,Ji_,Jf_,tranfreq_,igammai_,igammaf_,ileveli_,ilevelf_,acoef_,energyf,energyi,itermi_,itermf_) &
        !$omp&   schedule(static) shared(itermi,itermf)
        do iswap = 1,nselect
           !
           !
           tranfreq_= tranfreq(iswap)
           !
           !print*,"tranfreq_ = ",tranfreq_
           !
           !if ( tranfreq_<nu1 .or.  tranfreq_>nu2 ) cycle
           !
           Ji_ = Ji(iswap)
           Jf_ = Jf(iswap)
           !
           igammai_ = igammai(iswap)
           igammaf_ = igammaf(iswap)
           ileveli_ = ileveli(iswap)
           ilevelf_ = ilevelf(iswap)
           !
           acoef_ = acoef(iswap)
           !
           if (Jf_>jmax.or.Ji_>jmax) cycle
           !
           !if (igammaf_<1.or.igammai_<1) then
           !  write(6,"('negative igamma for iswap =',i8,'j,j,gam,gam=',4i6,'freq = ',f18.6)") jf_,igammaf_,ji_,igammai_,tranfreq_
           !endif 
           !
           energyf=level(jf_,igammaf_)%energy(ilevelf_)
           energyi=level(ji_,igammai_)%energy(ileveli_)
           !
           if (abs(energyf-energyi-tranfreq_)>1e-3) then
             write(6,"('energyf,energyi,nu,tranfreq= ',4f16.5)") energyf,energyi,energyf-energyi,tranfreq_
             write(6,"('jf,igammaf,ilevelf,ji,igammai,ileveli',6i9)") jf_,igammaf_,ilevelf_,ji_,igammai_,ileveli_
             stop 'wrong combination of states'
           endif
           !
           itermi_ = level(Ji_,igammai_)%nn(ileveli_)
           !
           itermf_ = level(Jf_,igammaf_)%nn(ilevelf_)
           !
           !print*,itermi_,itermf_,iswap
           !
           !
           !ilevelf_ = 0
           !!
           !do j = 1,nlevels
           !  !
           !  if (jf ==Jktau(j,1) .and. igammaf_ == sym(j) .and. abs(energyf-energies(j)<1e-4)) then 
           !    !
           !    ilevelf_ = j
           !    !
           !    exit
           !    !
           !  endif
           !  !
           !enddo
           !
           if (itermi_==0.or.itermf_==0) then 
             !
             print("( 'cannot find the correspondence in the energy-file for ',i,'th transition, levels ',i8,' and ',i8,', energies: ',2f14.4 )"), i,ileveli_,ilevelf_,energyi,energyf
             stop 'cannot find the correspondence in the energy-file'
             !
             itermf_ = -1 ; itermi_ = -1
             !
             cycle
             !
           endif
           !
           itermi(iswap) = itermi_
           itermf(iswap) = itermf_
           !
        enddo
        !$omp end parallel do
        !
        print*,"Printing to trans-file...",kswaps,nselect
        !
        write(unit_f,"(2i12,2x,es16.8,2x,f16.6)"),(itermf(iswap),itermi(iswap),acoef(iswap),tranfreq(iswap),iswap=1,nselect)
        !
        !!!write(unit_f,"(2i12,2x,es16.8)"),(itermf(iswap),itermi(iswap),acoef(iswap),iswap=1,nselect)
        !
        !call flash(unit_f)
        !
        !do iswap = 1, nselect
           !
           !tranfreq_= tranfreq(iswap)
           !
           !if ( tranfreq_<nu1 .or.  tranfreq_>nu2 ) cycle
           !
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !if (i<=117.and.(level(jf,igammaf_)%energy(ilevelf_)-level(ji,igammai)%energy(ileveli_))>10000.0000)  cycle
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !
           !if (energyf>12000.0d0) cycle
           !unit_f = 2000+int(energyf/1000)
           !
           !
           !unit_f = 1000+int(tranfreq(iswap)/100)
           !
           !write(fchar, '(i4)') unit_f
           !fname = trim(intfilenameout)//'_'//trim(adjustl(fchar))//'.trans'
           !
           !INQUIRE (FILE = fname, opened = ifopen)
           !
           !if (.not.ifopen)  open(unit=unit_f,file=fname,action='write',status='replace')
           !
           !unit_f = 13
           !
           !write(unit_f,"(2i12,2x,es16.8,2x,f16.6)"),itermf(iswap),itermi(iswap),acoef(iswap),tranfreq(iswap)
           !
         !enddo
         !
         ilines = ilines + nselect
         !
         if (nswap_<nswap) exit loop_file
         !
      enddo loop_file
      !
      close(unit_f)
      !
      !if (imis /= 0) print("( 'missassigned idescr = ',i)"), imis
      !
    enddo 
    !
    do i = 0,jmax
     print("(i4,2x,i14)"), i,Nlines(i)
    enddo
    !
    deallocate(sym,energies,Jktau)
    !
    print("( 'Done!')")
    !
    close(1)
    close(13)


end program linelist_trove
