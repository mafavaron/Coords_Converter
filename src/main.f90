program cd2ll

    use COORDS
    use ll_utm

    implicit none
    
    ! Locals
    integer             :: iRetCode
    character(len=256)  :: sInputFile
    character(len=256)  :: sOutputFile
    character(len=256)  :: sBuffer
    type(Geo_Coord)     :: tLatLon
    type(UTM_Coord)     :: tUTM
    real(8)             :: rEast
    real(8)             :: rNorth
    real(8)             :: rConc
    
    ! Get parameters
    if(COMMAND_ARGUMENT_COUNT() /= 3) then
        print *, "cd2ll - Program converting from ConcDecode-generated fields (UTM) to Lat/Lon"
        print *
        print *, "Usage:"
        print *
        print *, "   cd2ll <Input_File> <UTM_Fuse> <Output_File>"
        print *
        print *, "Copyright 2019 by Servizi Territorio srl"
        print *, "                  All rights reserved"
        print *
        print *, "Written by Mauri Favaron"
        print *
        stop
    end if
    call GET_COMMAND_ARGUMENT(1, sInputFile)
    call GET_COMMAND_ARGUMENT(2, tUTM % Zone)
    call GET_COMMAND_ARGUMENT(3, sOutputFile)
    
    ! Get input file
    open(10, file=sInputFile, status='old', action='read', iostat=iRetCode)
    if(iRetCode /= 0) then
        print *, 'cd2ll:: error: Input file not found'
        stop
    end if
    read(10, "(a)", iostat=iRetCode) sBuffer
    if(iRetCode /= 0) then
        print *, 'cd2ll:: error: Input file is empty'
        stop
    end if
    open(11, file=sOutputFile, status='unknown', action='write')
    write(11, "('X, Y, Conc')")
    do
        read(10, *, iostat=iRetCode) rEast, rNorth, rConc
        if(iRetCode /= 0) exit
        tUTM % E = rEast
        tUTM % N = rNorth
        call UTMtoLL(WGS_84, tUtm, tLatLon)
        write(11, "(2(f16.13,','),e15.7)") tLatLon % Lon, tLatLon % Lat, rConc
    end do
    close(11)
    close(10)

end program cd2ll
