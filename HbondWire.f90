program analyse

implicit none
integer :: i, j, k, f, Nsolute, Nwat, Natoms, Nframes=100000
integer, dimension(2) :: point
real, dimension(3) :: box, dX 
real,allocatable :: X(:,:)
real, dimension(3) :: rOO, rOiH, rOjH
real :: mod_rOO, mod_rOiH, mod_rOjH
real, dimension(:,:), allocatable :: Hb
character(len=2),allocatable :: atom(:)
character(len=80) :: trajfile, title, inputfile

call getarg(1,inputfile)
open(10,file=inputfile)
read(10,*) Nsolute
read(10,*) point
read(10,*) box
read(10,*) trajfile
open(11,file=trajfile)
Nwat=Nsolute+1

do f=1,Nframes

  read(11,*,end=100) Natoms
  read(11,'(A50)') title
  if(f.EQ.1) allocate(X(Natoms,3))
  if(f.EQ.1) allocate(atom(Natoms))
  if(f.EQ.1) allocate(Hb(Natoms,Natoms))
  do i=1,Natoms
    read(11,*) atom(i),X(i,:)
  enddo

!wrap water
  do i=Nwat,Natoms,3
    dX(:)=X(i,:)-X(point(1),:)
    X(i,:) = X(i,:) - box(:)*ANINT(dX(:)/box(:))
    X(i+1,:) = X(i+1,:) - box(:)*ANINT(dX(:)/box(:))
    X(i+2,:) = X(i+2,:) - box(:)*ANINT(dX(:)/box(:))
  enddo

!H-bonds:
!between water molecules
  Hb = 0.
  do i=Nwat,Natoms,3
    do j=i,Natoms,3
      if(j .NE. i) then
        rOO(:)=X(i,:)-X(j,:)
        mod_rOO = sqrt(dot_product(rOO,rOO))
        do k=Nwat,Natoms
          if(atom(k) .EQ. "H ") then
            rOiH(:)=X(i,:)-X(k,:)
            mod_rOiH = sqrt(dot_product(rOiH,rOiH))
            rOjH(:)=X(j,:)-X(k,:)
            mod_rOjH = sqrt(dot_product(rOjH,rOjH))
            if(mod_rOiH .lt. mod_rOjH) then
              Hb(i,j) = Hb(i,j) + (1 - ((mod_rOiH+mod_rOjH-mod_rOO)/0.6)**8) / (1 - ((mod_rOiH+mod_rOjH-mod_rOO)/0.6)**12)
            else
              Hb(j,i) = Hb(j,i) + (1 - ((mod_rOiH+mod_rOjH-mod_rOO)/0.6)**8) / (1 - ((mod_rOiH+mod_rOjH-mod_rOO)/0.6)**12)
            endif
          endif
        enddo
        Hb(i,j) = Hb(i,j) * (1 - ((mod_rOO - 2.7)/0.5)**10) / (1 - ((mod_rOO - 2.7)/0.5)**20)
        Hb(j,i) = Hb(j,i) * (1 - ((mod_rOO - 2.7)/0.5)**10) / (1 - ((mod_rOO - 2.7)/0.5)**20)
      endif
    enddo
  enddo

!between water molecules and end-point_1 (nucleophilic water)
  do j=Nwat,Natoms,3
    rOO(:)=X(point(1),:)-X(j,:)
    mod_rOO = sqrt(dot_product(rOO,rOO))
    do k=Nwat,Natoms
      if(atom(k) .EQ. "H ") then
        rOiH(:)=X(point(1),:)-X(k,:)
        mod_rOiH = sqrt(dot_product(rOiH,rOiH))
        rOjH(:)=X(j,:)-X(k,:)
        mod_rOjH = sqrt(dot_product(rOjH,rOjH))
        if(mod_rOiH.lt.1.3 .or. mod_rOjH.lt.1.3) then
          if(mod_rOiH .lt. mod_rOjH) then
            Hb(point(1),j) = Hb(point(1),j) + &
             (1 - ((mod_rOiH+mod_rOjH-mod_rOO)/0.6)**8) / (1 - ((mod_rOiH+mod_rOjH-mod_rOO)/0.6)**12)
          else
            Hb(j,point(1)) = Hb(j,point(1)) + &
             (1 - ((mod_rOiH+mod_rOjH-mod_rOO)/0.6)**8) / (1 - ((mod_rOiH+mod_rOjH-mod_rOO)/0.6)**12)
          endif
        endif
      endif
    enddo

    do k=point(1)+1, point(1)+2
      if(atom(k) .EQ. "H ") then
        rOiH(:)=X(point(1),:)-X(k,:)
        mod_rOiH = sqrt(dot_product(rOiH,rOiH))
        rOjH(:)=X(j,:)-X(k,:)
        mod_rOjH = sqrt(dot_product(rOjH,rOjH))
        if(mod_rOiH.lt.1.3 .or. mod_rOjH.lt.1.3) then
          if(mod_rOiH .lt. mod_rOjH) then
            Hb(point(1),j) = Hb(point(1),j) + &
             (1 - ((mod_rOiH+mod_rOjH-mod_rOO)/0.6)**8) / (1 - ((mod_rOiH+mod_rOjH-mod_rOO)/0.6)**12)
          else
            Hb(j,point(1)) = Hb(j,point(1)) + &
             (1 - ((mod_rOiH+mod_rOjH-mod_rOO)/0.6)**8) / (1 - ((mod_rOiH+mod_rOjH-mod_rOO)/0.6)**12)
          endif
        endif
      endif
    enddo

    Hb(point(1),j) = Hb(point(1),j) * &
         (1 - ((mod_rOO - 2.7)/0.5)**10) / (1 - ((mod_rOO - 2.7)/0.5)**20)
    Hb(j,point(1)) = Hb(j,point(1)) * &
         (1 - ((mod_rOO - 2.7)/0.5)**10) / (1 - ((mod_rOO - 2.7)/0.5)**20)
  enddo

!between water molecules and end-point_2
  do j=Nwat,Natoms,3
    rOO(:)=X(point(2),:)-X(j,:)
    mod_rOO = sqrt(dot_product(rOO,rOO))
    do k=Nwat,Natoms
      if(atom(k) .EQ. "H ") then
        rOiH(:)=X(point(2),:)-X(k,:)
        mod_rOiH = sqrt(dot_product(rOiH,rOiH))
        rOjH(:)=X(j,:)-X(k,:)
        mod_rOjH = sqrt(dot_product(rOjH,rOjH))
        if(mod_rOiH.lt.1.3 .or. mod_rOjH.lt.1.3) then
          if(mod_rOiH .lt. mod_rOjH) then
            Hb(point(2),j) = Hb(point(2),j) + &
             (1 - ((mod_rOiH+mod_rOjH-mod_rOO)/0.6)**8) / (1 - ((mod_rOiH+mod_rOjH-mod_rOO)/0.6)**12)
          else
            Hb(j,point(2)) = Hb(j,point(2)) + &
             (1 - ((mod_rOiH+mod_rOjH-mod_rOO)/0.6)**8) / (1 - ((mod_rOiH+mod_rOjH-mod_rOO)/0.6)**12)
          endif
        endif
      endif
    enddo
    Hb(point(2),j) = Hb(point(2),j) * &
         (1 - ((mod_rOO - 3.7)/0.5)**10) / (1 - ((mod_rOO - 3.7)/0.5)**20)
    Hb(j,point(2)) = Hb(j,point(2)) * &
         (1 - ((mod_rOO - 3.7)/0.5)**10) / (1 - ((mod_rOO - 3.7)/0.5)**20)
  enddo


  open(21,file="nuc_wat_chains ")
  if(Hb(point(1),point(2)) .GE. 0.8) then
    write(21,*) f, i, Hb(point(1),point(2))
  endif
  open(22,file="1_water_chains")
  open(23,file="2_water_chains")
  open(24,file="3_water_chains")
  do i=Nwat,Natoms,3
    if(Hb(point(1),i)*Hb(i,point(2)) .GE. 0.8) then
      write(22,*) f, i, Hb(point(1),i), Hb(i,point(2))
    endif
    do j=Nwat,Natoms,3
      if(Hb(point(1),i)*Hb(i,j)*Hb(j,point(2)) .GE. 0.8) then
        write(23,*) f, i, j, Hb(point(1),i), Hb(i,j), Hb(j,point(2))
      endif
      do k=Nwat,Natoms,3
        if(Hb(point(1),i)*Hb(i,j)*Hb(j,k)*Hb(k,point(2)) .GE. 0.8) then
         write(24,*) f, i, j, k, Hb(point(1),i), Hb(i,j), Hb(j,k), Hb(k,point(2))
        endif
      enddo !k
    enddo !j
  enddo !i

enddo !f
100 continue
 print*,"last frame",f-1

end program
