!!------------------------------------------------------------------------------
!! k-d tree for fast retrieval of observations.
!!
!! Travis Sluka
!! University of Maryland
!! tsluka@umd.edu
!!------------------------------------------------------------------------------




!!------------------------------------------------------------
!!------------------------------------------------------------
!! Kd-tree
!!------------------------------------------------------------
module KDTREE
  use common_obs
  implicit none
  private

  public :: kd_root, kd_init, kd_search


  !! a node for the kd-tree structure
  type boxnode
     real(r_size) :: lo(2), hi(2)  !! high and low of lon/lat
     integer :: mom, dau1, dau2, ptlo, pthi
  end type boxnode

  !! wrapper for all the information needed by the kd-tree after
  !! creation and during searching. User shouldn't have to
  !! access any of these members directly
  type KD_ROOT
    integer, pointer :: ptindx(:)
    type(boxnode), pointer :: boxes(:)
    real(r_size) :: lon_bounds(2)
  end type KD_ROOT

contains




  !!------------------------------------------------------------
  !! Initialize a Kd Tree structure.
  !! TODO, probably is not currently optimized for dealing with
  !!  observations with identical lat/lon. More boxes are being
  !!  generated than needed in that case.
  !!------------------------------------------------------------
  subroutine kd_init(root, lons, lats)
    type(KD_ROOT), intent(out) :: root       !! the root node of our kd-tree
    real(r_size), target :: lons(:), lats(:)

    !! variables for kd-tree creation loop
    integer :: kk, np, ntmp, n, m, nboxes, nowtask, ptlo, pthi
    integer :: tmom, tdim, jbox
    real(r_size), pointer :: cp(:)
    integer,pointer :: hp(:)
    real(r_size) :: lo(2), hi(2)
    integer :: taskmom(50), taskdim(50)


    !! generate initial unsorted index array
    write (*,*) "Initializting k-d tree with ",size(lons), "observations"
    allocate(root%ptindx(size(lons)))
    do n=1,size(root%ptindx)
       root%ptindx(n) = n
    end do

    !! save the longitude bounds, we'll need those in the future
    root%lon_bounds(1) = minval(lons)
    root%lon_bounds(2) = maxval(lons)
    if (root%lon_bounds(2)-root%lon_bounds(1) > 360) then
       write (*,*) "ERROR: kd_init(), the range of longitudes of the observations is too great (",&
            root%lon_bounds(1), root%lon_bounds(2),") kd_init doesn't care what the values are, just as long ",&
            " as they only span 360 degrees"
       stop 1
    end if

    !! calculate number of kd boxes needed and create memory for them
    m = 1
    ntmp = size(root%ptindx)
    do while (ntmp > 0)
       ntmp = ishft(ntmp,-1)
       m = ishft(m,1)
    end do
    nboxes = 2*size(root%ptindx)-ishft(m,-1)
    if (m<nboxes) nboxes = m
    nboxes = nboxes - 1
    allocate(root%boxes(nboxes))

    !! initialize the root box, and put its subdivision on the task list
    lo = (/-1e4, -1e4/)
    hi = (/1e4, 1e4/)
    root%boxes(1) = boxnode(lo,hi,0,0,0,1,size(root%ptindx))
    jbox = 1
    taskmom(1) = 1
    taskdim(1) = 0
    nowtask = 1

    do while(nowtask > 0)
       !! get the box sitting on the top of the task list
       tmom = taskmom(nowtask)
       tdim = taskdim(nowtask)
       nowtask = nowtask - 1
       ptlo = root%boxes(tmom)%ptlo
       pthi = root%boxes(tmom)%pthi
       hp => root%ptindx(ptlo:pthi)

       !! alternate between dividing by latitude and longitude
       if (tdim == 0) then
          cp => lons
       else
          cp => lats
       end if

       !! determine dividing points
       np = pthi - ptlo +1  !! total points
       kk = (np+1)/2        !! left most point of first subdivision

       !! do the array partitioning
       call kd_selecti(kk, hp, cp)

       !! create the daughters and push them onto the task
       !! list if they need further subdivision
       hi = root%boxes(tmom)%hi
       lo = root%boxes(tmom)%lo
       lo(tdim+1) = cp(hp(kk))
       hi(tdim+1) = lo(tdim+1)
       root%boxes(jbox+1)= boxnode(root%boxes(tmom)%lo,hi, tmom,0,0,ptlo, ptlo+kk-1)
       root%boxes(jbox+2)= boxnode(lo,root%boxes(tmom)%hi, tmom,0,0,ptlo+kk, pthi)
       jbox = jbox + 2
       root%boxes(tmom)%dau1 = jbox-1
       root%boxes(tmom)%dau2 = jbox

       !! subdivide left further
       if (kk > 2) then
          nowtask = nowtask +1
          taskmom(nowtask) = jbox-1
          taskdim(nowtask) = mod((tdim+1),2)
       end if
       !! subdivie right further
       if (np-kk > 4) then
          nowtask = nowtask + 1
          taskmom(nowtask) = jbox
          taskdim(nowtask) = mod((tdim+1),2)
       end if
    end do
  end subroutine kd_init





  subroutine kd_search(root, lons, lats, s_point, s_radius, r_points, r_distance, r_num)
    type(KD_ROOT), intent(in)    :: root
    real(r_size),  intent(in)    :: lats(:), lons(:)
    real(r_size),  intent(in)    :: s_point(2)         !! center of search circle
    real(r_size),  intent(in)    :: s_radius           !! radius, in m of search circle
    integer,       intent(inout) :: r_points(:)
    real(r_size),  intent(inout) :: r_distance(:)
    integer,       intent(out)   :: r_num

    real(r_size) :: r, dl, mid
    real(r_size) :: s_box_min(2), s_box_max(2), s_pt(2)
    integer :: r_num2, i
    
    !! some basic checks
    if (size(r_points) /= size(r_distance) ) then
       write (*,*) "ERROR: kd_search(), r_points and r_distance must be allocated with same size"
       stop 1
    end if

    !! determine the rough lat/lon bounding box for the search radius
    r = s_radius / re
    !! the min/max lat is easy...
    s_box_min(2) = s_point(2) - (r * 180.0 / pi)
    s_box_max(2) = s_point(2) + (r * 180.0 / pi)
    if (s_box_min(2) < -90) s_box_min(2) = -90
    if (s_box_max(2) >  90) s_box_max(2) = 90
    !! the min/max lon is more annoying
    if (s_box_max(2) == 90.0 .or. s_box_min(2) == -90.0) then
       s_box_min(1) = s_point(1) - 180.0
       s_box_max(1) = s_point(1) + 180.0
    else
       dl = asin(sin(r) / cos(s_point(2)*pi/180.0)) * 180.0 / pi
       s_box_min(1) = s_point(1) - dl
       s_box_max(1) = s_point(1) + dl
    end if

    !! find points given the supplied lat/lon
    r_num=0
    call kd_search_inner(root,lons,lats,s_point,s_radius,r_points,r_distance,r_num)

    !! check to see if our search area extends past the longitude of the given obs,
    !! if so we need to run the search again with the search area shifted by 360 degrees
    r_num2 = 0
    mid = (root%lon_bounds(1)+root%lon_bounds(2)) * 0.5 !! middle longitude of the observations
    s_pt = s_point
    if (s_point(1) < mid .and. s_box_min(1) < root%lon_bounds(1)) then
       !! search area overlaps past the left of the map
       s_pt(1) = s_pt(1)+360
       call kd_search_inner(root,lons,lats,s_pt,s_radius,r_points(r_num+1:),r_distance(r_num+1:),r_num2)
    else if (s_point(1) > mid .and. s_box_max(1) > root%lon_bounds(2)) then
       !! search area overlaps past the right of the map
       s_pt(1) = s_pt(1)-360
       call kd_search_inner(root,lons,lats,s_pt,s_radius,r_points(r_num+1:),r_distance(r_num+1:),r_num2)
    end if

    !! if we ran 2 search, and there is the possibility of duplicates (longitudes wrap around and touch)...
    r_num = r_num+r_num2
    write (*,*) s_box_min, s_box_max
    if (r_num2 > 0 )then !.and. (360+s_box_min(1) - s_box_max(1) < 1)) then
       !! sort the list and remove duplicates
       call qsort(r_points(1:r_num))
       do i=r_num-1,1,-1
          if( r_points(i) == r_points(i+1) ) then
             call swap(r_points(i), r_points(r_num))
             r_num = r_num-1
          end if
       end do
    end if
    
  end subroutine kd_search


  
  !! --------------------------------------------------------------------------------
  !! The inner loop of the kd search, does the heavy lifting of traversing the tree
  !! --------------------------------------------------------------------------------
  subroutine  kd_search_inner(root, lons, lats, s_point, s_radius, r_points, r_distance, r_num)
    type(KD_ROOT), intent(in)    :: root
    real(r_size),  intent(in)    :: lats(:), lons(:)
    real(r_size),  intent(in)    :: s_point(2)         !! center of search circle
    real(r_size),  intent(in)    :: s_radius           !! radius, in m of search circle
    integer,       intent(inout) :: r_points(:)
    real(r_size),  intent(inout) :: r_distance(:)
    integer,       intent(out)   :: r_num

    real(r_size) :: r, dl
    real(r_size) :: s_box_min(2), s_box_max(2)
    integer :: k, i, n, nb, nbold, ntask, jdim, d1, d2
    integer :: task(50)

    real(r_size) :: latr, slatr,clatr
    real(r_size) :: plat, plon
    type(boxnode), pointer :: boxes(:)

    boxes => root%boxes

    !! determine the rough lat/lon bounding box for the search radius
    r = s_radius / re
    !! the min/max lat is easy...
    s_box_min(2) = s_point(2) - (r * 180.0 / pi)
    s_box_max(2) = s_point(2) + (r * 180.0 / pi)
    if (s_box_min(2) < -90) s_box_min(2) = -90
    if (s_box_max(2) >  90) s_box_max(2) = 90
    !! the min/max lon is more annoying
    if (s_box_max(2) == 90.0 .or. s_box_min(2) == -90.0) then
       s_box_min(1) = s_point(1) - 180.0
       s_box_max(1) = s_point(1) + 180.0
    else
       dl = asin(sin(r) / cos(s_point(2)*pi/180.0)) * 180.0 / pi
       s_box_min(1) = s_point(1) - dl
       s_box_max(1) = s_point(1) + dl
    end if
    

    !! find the smallest box that completely contains the bounds of the search point
    nb = 1
    jdim = 0
    do while(boxes(nb)%dau1 /= 0)
       nbold = nb
       d1 = boxes(nb)%dau1
       d2 = boxes(nb)%dau2
       if ( s_box_max(jdim+1) <= boxes(d1)%hi(jdim+1) ) then
          nb = d1
       else if (  s_box_min(jdim+1) >= boxes(d2)%lo(jdim+1) ) then
          nb = d2
       end if
       jdim = mod(jdim+1, 2)
       if (nb == nbold) exit
    end do

    !! traverse the tree below this starting box, only as needed
    task(1) = nb
    ntask = 1
    r_num = 0
    clatr = cos(s_point(2)*pi/180)
    slatr = sin(s_point(2)*pi/180)    

    plat = -100
    plon = -1000

    do while(ntask /= 0)
       k = task(ntask)
       ntask = ntask - 1

       !! ignore boxes definitely outside the radius
       i = 0
       do n = 1,2
          if( (boxes(k)%lo(n) - s_box_max(n)) * (boxes(k)%hi(n)-s_box_min(n)) > 0 ) then
             i = 1
             exit
          end if
       end do
       if (i == 1) cycle

       
       if(boxes(k)%dau1 /= 0) then
          !! process child boxes
          task(ntask+1) = boxes(k)%dau1
          task(ntask+2) = boxes(k)%dau2
          ntask = ntask + 2
       else
          !! process points in this box
          do i=boxes(k)%ptlo, boxes(k)%pthi
             n = root%ptindx(i)
             !! don't bother calculating the radius if the previous point was the same.
             !! This case will happen frequently as there are often profiles at a same lat/lon
             if ( abs(plat-lats(n)) < 1.0d-4 .and. abs(plon-lons(n)) < 1.0d-4) then
             else
                plat = lats(n)
                plon=lons(n)
                latr = lats(n)*pi/180
                ! TODO, this line is very expensive, can we calulate distance with less accurate but faster trig functions?
                r = re*acos(slatr*sin(latr) + clatr*cos(latr)*cos( (lons(n)-s_point(1))*pi/180.0 ))
             end if
             
             if (r < s_radius) then
                !! make sure there's room for new found obs
                if (r_num == size(r_points)) then
                   write (*,*) "ERROR: more observations found than there was room for in kd_search(). Increase size of the ",&
                        " r_points and r_distance arrays"
                   stop 1
                end if
                !! add the obs
                r_num = r_num + 1
                r_points(r_num) = n
                r_distance(r_num) = r
             end if             
          end do
       end if
    end do
  end subroutine kd_search_inner




  !! --------------------------------------------------------------------------------
  !! permutes indx[1...n] to make
  !!   arr[indx[1..k-1]] <= arr[indx[k]] <= arr[indx[k+1,n]]
  !!   the array "arr" is not modified
  !! --------------------------------------------------------------------------------
  subroutine kd_selecti(k, indx, arr)
    integer, intent(in) :: k
    integer, intent(inout) :: indx(:)
    real(r_size),intent(in) :: arr(:)

    integer :: i,ia,ir,j,l,mid
    real(r_size) :: a

    l=1
    ir=size(indx)
    do while(.true.)
       if (ir <= l+1) then
          if (ir == l+1 .and. arr(indx(ir)) < arr(indx(l))) &
             call swap(indx(l), indx(ir))
          exit
       else
          mid = (l+ir) / 2
          call swap(indx(mid), indx(l+1))
          if (arr(indx(l)) > arr(indx(ir)))   call swap(indx(l),indx(ir))
          if (arr(indx(l+1)) > arr(indx(ir))) call swap(indx(l+1),indx(ir))
          if (arr(indx(l)) > arr(indx(l+1)))  call swap(indx(l),indx(l+1))
          i=l+1
          j=ir
          ia=indx(l+1)
          a=arr(ia)
          do while(.true.)
             i = i +1
             do while(arr(indx(i)) < a)
                i = i + 1
             end do
             j = j-1
             do while(arr(indx(j)) > a)
                j = j - 1
             end do
             if (j < i) exit
             call swap(indx(i), indx(j))
          end do
          indx(l+1)=indx(j)
          indx(j)=ia
          if (j >= k) ir=j-1
          if (j<= k) l = i
       end if
    end do
  end subroutine kd_selecti




  !! --------------------------------------------------------------------------------
  !! Quicksort
  !!  Performs a standard quicksort, with fallback to insertion sort for the smallest
  !!  array segments. The array "arr" is an array of indices that are sorted using "key"
  !!  as the sorting key value. Pivot points are chosen by a median of 3 using the
  !!  segment start, middle, and end... if anyone cares.
  !! --------------------------------------------------------------------------------
  subroutine qsort(arr)
    integer, intent(inout) :: arr(:)        !! the array of indices to key to sort

    integer, parameter :: M=7          !! the max array size before switching over to insertion sort
    integer, parameter :: NSTACK=64    !! size of stack, should be big enough, increase if problems
    integer :: n, ir, a, l, i, j, k, jstack
    integer :: istack(NSTACK)

    !! initialize the pointers
    n = size(arr)
    ir = n
    l = 1
    jstack=-1

    do while (.true.)
       !! if the array segment is small enough, do an insertion sort
       !! ------------------------------------------------------------
       if (ir-l < M) then
          !! quick and dirty insertion sort
          do j=l+1,ir
             a = arr(j)
             do i=j-1,l,-1
                if (arr(i) < a) exit
                arr(i+1) = arr(i)
             end do
             arr(i+1) = a
          end do
          !! if there are no other segments left to sort, then we are done!
          if (jstack < 0) exit
          !! pop the stack and get the next segment for partitioning
          ir=istack(jstack)
          l=istack(jstack-1)
          jstack = jstack- 2

       !! Otherwise, for a large segment, do the quicksort
       !!------------------------------------------------------------
       else
          !! pick the left, center, and right elements of the segment,
          !! find the median as the partion point. make sure those 3
          !! elements are in order
          k= ishft((l+ir),-1)
          call swap(arr(k), arr(l+1))
          if (arr(l)   > arr(ir))  call swap(arr(l), arr(ir))
          if (arr(l+1) > arr(ir))  call swap(arr(l+1), arr(ir))
          if (arr(l)   > arr(l+1)) call swap(arr(l), arr(l+1))
          !! initialize pointers for partitioning
          i = l+1
          j = ir
          a = arr(l+1) !! index of partitioning element
          do while(.true.)
             !! scan up to find element > a
             i = i +1
             do while(arr(i) < a)
                i = i + 1
             end do
             !! scan down to find element < a
             j = j-1
             do while(arr(j) > a)
                j = j - 1
             end do
             !! pointers have crossed, partitioning is done...
             if (j < i) exit
             !! otherwise, swap elements and continue partitioning
             call swap(arr(i), arr(j))
          end do
          !! put the partitioning elment back in
          arr(l+1) = arr(j)
          arr(j)=a
          !! push pointers of larger array onto stack, process smaller
          !! array immediately
          jstack = jstack + 2
          if (jstack >= NSTACK) then
             write(*,*) "NSTACK is too small!!"
             stop 1
          end if
          if (ir-i+1 >= j-1) then
             istack(jstack) = ir
             istack(jstack-1) = i
             ir=j-1
          else
             istack(jstack) = j-1
             istack(jstack-1) = l
             l=i
          end if
       end if
    enddo
  end subroutine qsort


  !! Convenience function to swap two array indices
  subroutine swap(a1, a2)
    integer, intent(inout) :: a1, a2
    integer :: a
    a = a1
    a1 = a2
    a2 = a
  end subroutine swap

end module KDTREE
