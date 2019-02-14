program pr
!-------------------------------------------------------------------------------
! f90
!-------------------------------------------------------------------------------
! moduli:
 use coords
 use GAUSS_BOAGA
 use ll_utm
 use LatLon
!-------------------------------------------------------------------------------
 implicit none

!-------------------------------------------------------------------------------
! coordinate conversion
 type(GEO_COORD)   :: tLL   ! tLL%Lat  tLL%Lon        (r8 r8)
 type(GAUSS_COORD) :: tGB   ! tGB%S  tGB%N  tGB%E     (i r8 r8)
 type(UTM_COORD)   :: tUTM  ! tUTM%Zone tUTM%N tUTM%E (c4 r8 r8)

!-------------------------------------------------------------------------------
 character(1) ans
 integer ians
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
 write (6,*) 'scegli input:'
 write (6,*) ' 1 : LonLat --> GB, UTM'
 write (6,*) ' 2 :     GB --> LonLat, UTM'
 write (6,*) ' 3 :    UTM --> LonLat, GB'
 read (5,'(i1)') ians

!-------------------------------------------------------------------------------
 if (ians == 1) then  ! ' 1 : LonLat --> GB, UTM'
  write (6,*) 'Longitude=?'
  read (5,*) tLL%Lon
  write (6,*) 'Latitude=?'
  read (5,*) tLL%Lat
  write (6,'(2(A,F20.8))') ' Lon=', tLL%Lon, '    Lat=', tLL%Lat

  tGB = LL_TO_GB(tLL)
  write (6,'(2(A,F20.8),A,I4)') ' GBx=', tGB%E, '    GBy=', tGB%N &
&    , ' tGB%S=WEST=', tGB%S

  call LLtoUTM(WGS_84, tLL, tUTM)
  write (6,'(2(A,F20.8),3A)') 'UTMx=',tUTM%E, '   UTMy=',tUTM%N &
&                         , ' UTM Zone: ''',tUTM%Zone,''''

 end if

!-------------------------------------------------------------------------------
 if (ians == 2) then  ! ' 2 :     GB --> LonLat, UTM'
  write (6,*) 'GBx=?'
  read (5,*) tGB%E
  write (6,*) 'GBy=?'
  read (5,*) tGB%N
  tGB%S= WEST
  write (6,'(2(A,F20.8),A,I4)') ' GBx=', tGB%E, '    GBy=', tGB%N &
&    , ' tGB%S=WEST=', tGB%S

  tLL = GB_TO_LL(tGB)
  write (6,'(2(A,F20.8))') ' Lon=', tLL%Lon, '    Lat=', tLL%Lat

  tUTM%Zone= '+32T'
 !write (6,*) 'tUTM%Zone=[', tUTM%Zone,']'

  call LLtoUTM(WGS_84, tLL, tUTM)
  write (6,'(2(A,F20.8),3A)') 'UTMx=',tUTM%E, '   UTMy=',tUTM%N &
&                         , ' UTM Zone: ''',tUTM%Zone,''''

 end if

!-------------------------------------------------------------------------------
 if (ians == 3) then  ! ' 3 :    UTM --> LonLat, GB'
  write (6,*) 'UTMx=?'
  read (5,*) tUTM%E
  write (6,*) 'UTMy=?'
  read (5,*) tUTM%N
  tUTM%Zone= '+32T'
  write (6,'(2(A,F20.8),3A)') 'UTMx=',tUTM%E, '   UTMy=',tUTM%N &
&                         , ' UTM Zone: ''',tUTM%Zone,''''

  call UTMtoLL(WGS_84, tUTM, tLL)
  write (6,'(2(A,F20.8))') ' Lon=', tLL%Lon, '    Lat=', tLL%Lat

  tGB = LL_TO_GB(tLL)
  write (6,'(2(A,F20.8),A,I4)') ' GBx=', tGB%E, '    GBy=', tGB%N &
&    , ' tGB%S=WEST=', tGB%S

 end if


!-------------------------------------------------------------------------------
 stop
end program pr
!-------------------------------------------------------------------------------
