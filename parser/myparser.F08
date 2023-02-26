
! A RUDIMENTARY MODULE TO PARSE COMMAND LINE ARGUMENTS
! CURRENT VERSION ACCEPTS ONLY A SET NUMBER OF FLAGS (NO OPTIONAL FLAGS)
! FLAGS CAN BE PRESENT IN ANY ORDER

! DESCRIPTION OF FLAGS
! -f <input file name>
! -o <log file name>

! STARTED ON 24 JUNE 2019
! UPDATED ON 05 JULY 2019

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! THE API !!!!!
    ! CHARACTER(len=len_filename), DIMENSION(maxflags) :: strinp
    ! CALL INPUT_FLAGS(strinp)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE my_cl_parser
    USE myglobal
    IMPLICIT NONE

    CONTAINS
    
    SUBROUTINE INPUT_FLAGS(strinp)
        CHARACTER(len=len_filename), DIMENSION(maxflags), INTENT(OUT)   :: strinp
        CHARACTER(len=len_filename), DIMENSION(maxflags*2) :: strarg
        CHARACTER(len=len_filename)                        :: tempstr
        INTEGER                                            :: myflags, i

        myflags = COMMAND_ARGUMENT_COUNT()
        IF(myflags .NE. maxflags*2) THEN
            WRITE(*,*) "ERROR) NEED ", maxflags, " ARGUMENTS (2*NUMFLAGS)"
            WRITE(*,*) "ERROR) RECEIVED ", maxflags, " ARGUMENTS"
            STOP
        END IF

        DO i = 1, maxflags*2
            CALL GET_COMMAND_ARGUMENT(i,strarg(i))
        END DO

        DO i = 1, maxflags
            tempstr = ADJUSTL(strarg(2*i-1))
            SELECT CASE(tempstr(1:2))
                CASE('-f','-F')
                    strinp(1) = ADJUSTL(strarg(2*i))
                CASE('-o','-O')
                    strinp(2) = ADJUSTL(strarg(2*i))
                CASE DEFAULT
                    WRITE(*,*) "ERROR) UNKNOWN FLAG "//tempstr(1:2)//" .. EXITING.."
                    STOP
            END SELECT
        END DO
    END SUBROUTINE INPUT_FLAGS
END MODULE my_cl_parser