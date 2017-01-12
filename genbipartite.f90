program genbipartite

! J.T. Kindt, Emory University, based on an algorithm developed in collaboration with C. Weeden

! Generates integer partitions of bipartite numbers (NA,NB), motivated by statistics of
! simulations of molecular clusters with two components.  

! These are represented as matrices with row and column indices ranging from 0
! to NA and 0 to NB; each matrix element represents the number of clusters whose
! composition is given by its row and column.  The sum of the product of all
! matrix elements with their row index will always equal NA, and likewise for
! the column index and NB.  

! Algorithm works by starting with the integer
! partitions of Nmax=NA+NB and working through all the ways of distributing A
! molecules among the distribution of cluster sizes corresponding to that
! integer partition.  

! Output is not in a regular "alphabetical", is grouped by the "parent" integer
! partitions of Nmax.


implicit none

integer, parameter:: Nmaxmax=80
! number of compositions seems to scale as ~0.5*exp(Nmax**0.66)
integer, parameter:: Ncompmax=100000

integer:: comp(Ncompmax,Nmaxmax)
integer:: nsum,ncompNmax,nbipartite
integer:: nocc(0:Nmaxmax,3),nlevels,nAocc(0:Nmaxmax)
integer:: current(0:Nmaxmax,0:Nmaxmax)
integer:: levelvec(3),nAlevel
integer:: nArem,inext,ilevel,n_nAk,mk
integer:: i,j,k,l,nc,nA,ndists,nsub,ipart
integer:: fill_limit,nmax,inextlevel,minocc,nrest,ncl
logical:: if_filled(0:Nmaxmax)

write (*,*) "number of monomers (maximum is ",Nmaxmax,")"
read (*,*) Nmax
write (*,*) "number of A + B = ",Nmax,". Enter number of A. "
read (*,*) nA
if (nA.gt.nMax) then
        write (*,*) "nA must be less than ",Nmax
        stop
endif

call genunipart (comp,ncompNmax,Nmax,Nmaxmax,Ncompmax)
!produces array comp(i,j) containing the number of j-mers in partition index i.

nbipartite=0
n_nAk=0
do ipart=1,ncompNmax  ! Loop over all integer partitions of Nmax
! reassign to levels
!    (this allows us to skip over values of j with no contribution to the
!    partition)

  write (*,*) "################# parent partition #",ipart,":", comp(ipart,1:Nmax) 
   nocc=0
   nlevels=0
   do j=1,Nmax
      if (comp(ipart,j).gt.0) then
              nlevels=nlevels+1
              nocc(nlevels,1)=j        ! keeps track of cluster size at next occupied level
              nocc(nlevels,2)=comp(ipart,j)!  occupation number of that level
              nocc(nlevels,3)=j*comp(ipart,j)  ! total # of monomers assoc. with that level
      endif
   enddo
!  fill from top
   nAocc=0
   inext=nlevels
   nArem=nA           ! nArem is # of A remaining, i.e. that have not been assigned;
                      !  initialized as all of them 
   if_filled=.FALSE.  ! if_filled is a logical array that is .TRUE. when
                      ! a level is entirely composed of A
   if_filled(0)=.TRUE.

   do while (inext.le.nlevels)


   do ilevel=inext,1,-1 ! distributing A molecules among the levels; fill from inext down
      nAocc(ilevel)=min(nArem,nocc(ilevel,3))  ! either add all of remaining A molecules to this level
                                               ! or fill it to its capacity
      nArem=nArem-nAocc(ilevel) 
      ! check if level is full
      if_filled(ilevel)=(nocc(ilevel,3).eq.nAocc(ilevel))  
   enddo

   ! start with max # of pure A k-mer clusters in each level; maximally
   ! segregated state
   current=0  ! current will contain the current partition in matrix form
   do ilevel=nlevels,1,-1
      nAlevel=nAocc(ilevel) 
      k=nocc(ilevel,1)
      mk=nocc(ilevel,2)
      call initlevel(nAlevel,k,mk,current,Nmaxmax)
   enddo
   
   nbipartite=nbipartite+1
   call outcurrent (current,nbipartite,Nmax,nA,Nmaxmax)  ! print the first partition in this series
       
   ! cycle through levels 

   ilevel=nlevels
   do while (ilevel.gt.0)

   do ilevel=nlevels,1,-1
      inextlevel=0
      nAlevel=nAocc(ilevel)
      k=nocc(ilevel,1)
      mk=nocc(ilevel,2)
      ! find minocc (minimum occupied, i.e. the cluster with total monomers k that has
      !           contains the fewest A molecules)
      do i=0,k
         if (current(i,k-i).gt.0) then
                 minocc=i
                 exit
         endif
      enddo
      inextlevel=k+1
      do i=minocc+2,k
         if (current(i,k-i).gt.0) then
                 inextlevel=i
                 exit
         endif
      enddo
      if (inextlevel.le.k) then ! can make a new distribution by taking an A
                                ! away from a cluster with #A=inextlevel and
                                ! redistributing A molecules from clusters with
                                ! fewer A's
          nrest=inextlevel
          ncl=1   ! ncl keeps count of the number of clusters
          current(inextlevel,k-inextlevel)=current(inextlevel,k-inextlevel)-1
          ! count all lower A's to distribute
          do j=inextlevel-1,0,-1
               nrest=nrest+j*current(j,k-j)
               ncl=ncl+current(j,k-j)
               current(j,k-j)=0
          enddo
          current(inextlevel-1,k-inextlevel+1)=nrest/(inextlevel-1)  ! load them
                                                                     ! into inextlevel-1
          ncl=ncl-current(inextlevel-1,k-inextlevel+1)
          current(mod(nrest,inextlevel-1),k-mod(nrest,inextlevel-1))=1  ! remainder get grouped together
          ncl=ncl-1
          current(0,k)=current(0,k)+ncl    ! any remaining clusters will have no A
          nbipartite=nbipartite+1
          call outcurrent(current,nbipartite,Nmax,nA,Nmaxmax)
          exit ! jump back and start over from highest level
      else
              ! if A molecules are distributed as evenly as possible
              ! among the clusters on this level, then reset this level and moving to the next
          call initlevel(nAlevel,k,mk,current,Nmaxmax)
      endif
   enddo

   enddo

   ! check lowest filled levels
   nArem=0
   fill_limit=0
   do ilevel=1,nlevels
      if (if_filled(ilevel)) then
            nArem=nArem+nAocc(ilevel)
            fill_limit=ilevel
      else
         exit
      endif
   enddo
   !write (*,*) "* fill_limit",fill_limit
   nArem=nArem+nAocc(fill_limit+1)  ! nArem contains the # of A's that have 
                                    ! all been distributed in filled lowest-k
                                    ! levels or the next highest level; i.e.,
                                    ! the # of A's that cannot still be moved down
                                    ! a level or more.
   !write (*,*) "* nArem",nArem 
   inext=nlevels+1
   do ilevel=fill_limit+2,nlevels
      if (nAocc(ilevel).gt.0) then
           nAocc(ilevel)=nAocc(ilevel)-1  ! take one A away from the next
                                          ! highest occupied level, start over
                                          ! distributing the nArem among the
                                          ! lower levels.
           inext=ilevel-1
           if_filled(ilevel)=.FALSE.
           exit
      endif
   enddo
   nArem=nArem+1
   !write (*,*) "nArem, inext at end of loop:",nArem,inext

   enddo ! while (inext.le.nlevels); if inext=nlevels+1, that means that the A
         !  molecules have been distributed to fill the smallest clusters, and
         !  it is time to move to the next unipartite partition.


enddo

end program 

!**************
  
subroutine initlevel (nA,k,mk,current,Nmax)

! distributes nA A molecules among mk k-mers with maximum segregation

integer,intent(in):: nA,mk,k,Nmax
integer,intent(inout):: current (0:Nmax,0:Nmax)
integer:: j


! clear any current values
do j=0,k
  current(j,k-j)=0
enddo
if (nA.eq.0) then      ! if there are no A molecules, all clusters are pure B
       current(0,k)=mk
else
       current(k,0)=nA/k  
       current (mod(nA,k),k-mod(nA,k))=1
       current(0,k)=current(0,k)+(mk-current(k,0)-1)
endif

return
end subroutine


subroutine outcurrent(current,nbipartite,Nmax,nA,Nmaxmax)

integer, intent(in):: Nmax,nA,Nmaxmax,current(0:Nmaxmax,0:Nmaxmax),nbipartite
integer:: i

write (*,*) "######",nbipartite
do i=0,nA
  write (*,*) current(i,0:Nmax-nA)
enddo

end subroutine

subroutine genunipart(comp,ncompNmax,Nmax,Nmaxmax,Ncompmax)

! inductive algorithm to generate partitions
! of integer i from set of partitions of i-1;
! simple but not particularly efficient.

integer, intent(in):: Nmax,Nmaxmax
integer, intent(out):: comp(Ncompmax,Nmaxmax),ncompNmax
integer:: complast(Ncompmax,Nmaxmax),ncomp(Nmaxmax)
integer:: i,j,k,l
i=1
comp=0
complast=0
complast(1,1)=1
ncomp(1)=1

!  Every partition of i that has at least one free monomer
!  can be generated from a partition of i-1 by adding a free monomer.
!  Every partition of i that has no free monomers, and whose smallest part is k+1,
!  can be generated from a partition of i-1 that has exactly 1 part of k and
!  adding 1 to it to make a part of size k+1.  
do i=2,Nmax
! write (*,*) "#ncomp(",i-1,")=",ncomp(i-1),". Starting ",i," of ",Nmax,"."
 comp=complast
 do j=1,ncomp(i-1)
    comp(j,1)=complast(j,1)+1  ! first add a monomer to each
 enddo
 ncomp(i)=ncomp(i-1)
 do j=1,ncomp(i-1)  ! loop over compositions of i-1
    do k=1,i-1 ! find first non-zero entry
       if (complast(j,k).gt.0) then 
          if (complast(j,k).ne.1) exit
          ! if the entry for the k-mer is 1, then 
          ! make a new composition by taking
          ! away a k-mer and adding a k+1-mer  
          ncomp(i)=ncomp(i)+1 ! increment # of compositions
          if (ncomp(i).gt.Ncompmax) then
                  write (*,*) i, ncomp(i),Ncompmax,"too many comps"
                  stop
          endif
          comp(ncomp(i),k+1)=complast(j,k+1)+1
          do l=k+2,i  ! remainder of composition stays the same
            comp(ncomp(i),l)=complast(j,l)
          enddo
          exit
       endif
    enddo
 enddo
 !write (*,*) i,ncomp(i)
 complast=comp
enddo

ncompNmax=ncomp(Nmax)
return

end subroutine
