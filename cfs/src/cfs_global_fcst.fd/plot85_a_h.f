         subroutine plot_85_a_sn(fld,fstr,lons,lats,levs,kdt,unitno)
           implicit none
           real fld(lons,levs,lats)
           real gg_ncar(lons,lats)
           integer lons,lats,levs,i,j,k,kdt,unitno
           character*(*) fstr
           real scale
  
          scale = 1.0

          write(unitno+kdt,1884) lons,lats,levs,0.0,0,0,0,0


          do k=1,levs  ,3
             do j=1,lats
                do i=1,lons
                   gg_ncar(i,j) = fld(i,k,lats-j+1)
                enddo
             enddo
           write(unitno+kdt,1881) k
           write(unitno+kdt,1882) fstr
           write(unitno+kdt,1885) gg_ncar*scale
          print*,'plot_85_a wrote map number k=',k
         enddo
c

1881      format (' sigma level ',i2)
1882      format (a)
1884      format (3i5,f9.3,4i5)
1885      format(6(e20.13))

        close(unitno+kdt)
        return
        end
c
         subroutine plot_85_a(fld,fstr,lons,lats,levs,kdt,unitno)
           implicit none
           real fld(lons,levs,lats)
           real gg_ncar(lons,lats)
           integer lons,lats,levs,i,j,k,kdt,unitno
           character*(*) fstr
           real scale
  
          scale = 1.0

          write(unitno+kdt,1884) lons,lats,levs,0.0,0,0,0,0


          do k=1,levs  ,3
             do j=1,lats
                do i=1,lons
                   gg_ncar(i,j) = fld(i,k,j)
                enddo
             enddo
           write(unitno+kdt,1881) k
           write(unitno+kdt,1882) fstr
           write(unitno+kdt,1885) gg_ncar*scale
          print*,'plot_85_a wrote map number k=',k
         enddo
c

1881      format (' sigma level ',i2)
1882      format (a)
1884      format (3i5,f9.3,4i5)
1885      format(6(e20.13))

        close(unitno+kdt)
        return
        end
c
         subroutine plot_85_h(fld,fstr,lons,lats,levs,kdt,unitno,yhalo)
           implicit none
           integer lons,lats,levs,i,j,k,kdt,unitno,yhalo
           real fld(lons,levs,lats)
           real gg_ncar(lons-3,lats-2*yhalo)
           integer lonl,latl
           character*(*) fstr
           real scale
  
          scale = 1.0
          lonl=lons-3
          latl=lats-2*yhalo

          write(unitno+kdt,1884) lonl,latl,levs,0.0,0,0,0,0


          do k=1,levs  ,3
             do j=1,latl
                do i=1,lonl
                   gg_ncar(i,j) = fld(1+i, levs+1-k ,yhalo+j)
                enddo
             enddo
           write(unitno+kdt,1881) k
           write(unitno+kdt,1882) fstr
           write(unitno+kdt,1885) gg_ncar*scale
          print*,'plot_85_h wrote map number k=',k
         enddo
c

1881      format (' sigma level ',i2)
1882      format (a)
1884      format (3i5,f9.3,4i5)
1885      format(6(e20.13))

        close(unitno+kdt)
        return
        end
         subroutine plot_85_exp(fld,fstr,lons,lats,levs,kdt,unitno)
           implicit none
           real fld(lons,levs,lats)
           real gg_ncar(lons,lats)
           integer lons,lats,levs,i,j,k,kdt,unitno
           character*(*) fstr
           real scale
                                                                                                       
          scale = 10.0
                                                                                                       
          write(unitno+kdt,1884) lons,lats,levs,0.0,0,0,0,0
                                                                                                       
                                                                                                       
          do k=1,levs  ,3
             do j=1,lats
                do i=1,lons
                   gg_ncar(i,j) =exp( fld(i,k,lats+1-j) )
                enddo
             enddo
           write(unitno+kdt,1881) k
           write(unitno+kdt,1882) fstr
           write(unitno+kdt,1885) gg_ncar*scale
          print*,'plot_85_a wrote map number k=',k
         enddo
c
                                                                                                       
1881      format (' sigma level ',i2)
1882      format (a)
1884      format (3i5,f9.3,4i5)
1885      format(6(e20.13))
                                                                                                       
        close(unitno+kdt)
        return
        end

