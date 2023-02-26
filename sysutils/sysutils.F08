
! A MODULE TO FACILITATE SYSTEM COMMANDS LIKE COPYING AND MOVING FILES AROUND

! STARTED ON 24 JUNE 2019
! UPDATED ON 05 JULY 2019

    ! INTEGER, PARAMETER :: len_parname = 10
    ! INTEGER, PARAMETER :: len_partype = 1
    ! INTEGER, PARAMETER :: len_filename = 400
    ! INTEGER, PARAMETER :: kind_real   = 8
    ! INTEGER, PARAMETER :: max_flags   = 1

! A SAFETY FEATURE CAN BE TO INDEX() FOR THE SUBSTRING '<HEADDIR>' WITHIN THE DIRECTORY NAMES PASSED TO ANY OF THE SUBROUTINES
! THIS WILL MAKE SURE THAT ALL THE OPERATIONS ARE OCCURING WITHIN THE HEADDIR AND NOT USING ANY FILES FROM OUTSIDE


MODULE sysutils
USE myglobal
IMPLICIT NONE

    CONTAINS

    SUBROUTINE FCOPY(fsource,ftarget)
        CHARACTER(LEN=len_filename), INTENT(IN) :: fsource
        CHARACTER(LEN=len_filename), INTENT(IN) :: ftarget
        CHARACTER(LEN=len_syscommand) :: mycommand, CMSG
        INTEGER :: EXSTAT, CSTAT, ICHK

        ICHK = INDEX(TRIM(ftarget),TRIM(HEADDIR))
        IF(ICHK .EQ. 0) THEN
            WRITE(*,*) "ERROR) FCOPY REQUESTED FOR PATH OUTSIDE HEADDIR"
            WRITE(*,*) "ERROR) HEADDIR : ", TRIM(HEADDIR)
            STOP
        END IF

        mycommand = 'cp -p '//TRIM(ADJUSTL(fsource))//' '//TRIM(ADJUSTL(ftarget))

        CALL EXECUTE_COMMAND_LINE(TRIM(mycommand),EXITSTAT=EXSTAT,CMDSTAT=CSTAT,CMDMSG=CMSG)

        IF(CSTAT .NE. 0) THEN
            WRITE(*,*) "ERROR) COULD NOT EXECUTE FCOPY COMMAND.."
            WRITE(*,*) "ERROR) COMMAND  : ", TRIM(mycommand)
            WRITE(*,*) "ERROR) EXIT MSG : ", TRIM(CMSG)
            STOP
        END IF

        IF(EXSTAT .NE. 0) THEN
            WRITE(*,*) "ERROR) EXECUTED FCOPY COMMAND THROWS AN ERROR.."
            WRITE(*,*) "ERROR) COMMAND  : ", TRIM(mycommand)
            STOP
        END IF
    END SUBROUTINE FCOPY

    SUBROUTINE DCOPY(dsource,dtarget)
        CHARACTER(LEN=len_filename), INTENT(IN) :: dsource
        CHARACTER(LEN=len_filename), INTENT(IN) :: dtarget
        CHARACTER(LEN=len_syscommand) :: mycommand, CMSG
        INTEGER :: EXSTAT, CSTAT, ICHK

        ICHK = INDEX(TRIM(dtarget),TRIM(HEADDIR))
        IF(ICHK .EQ. 0) THEN
            WRITE(*,*) "ERROR) DCOPY REQUESTED FOR PATH OUTSIDE HEADDIR"
            WRITE(*,*) "ERROR) HEADDIR : ", TRIM(HEADDIR)
            STOP
        END IF

        mycommand = 'cp -rp '//TRIM(ADJUSTL(dsource))//' '//TRIM(ADJUSTL(dtarget))

        CALL EXECUTE_COMMAND_LINE(TRIM(mycommand),EXITSTAT=EXSTAT,CMDSTAT=CSTAT,CMDMSG=CMSG)

        IF(CSTAT .NE. 0) THEN
            WRITE(*,*) "ERROR) COULD NOT EXECUTE DCOPY COMMAND.."
            WRITE(*,*) "ERROR) COMMAND  : ", TRIM(mycommand)
            WRITE(*,*) "ERROR) EXIT MSG : ", TRIM(CMSG)
            STOP
        END IF

        IF(EXSTAT .NE. 0) THEN
            WRITE(*,*) "ERROR) EXECUTED DCOPY COMMAND THROWS AN ERROR.."
            WRITE(*,*) "ERROR) COMMAND  : ", TRIM(mycommand)
            WRITE(*,*) "ERROR) EXIT MSG : ", TRIM(CMSG)
            STOP
        END IF
    END SUBROUTINE DCOPY

    SUBROUTINE FMOVE(fsource,ftarget)
        CHARACTER(LEN=len_filename), INTENT(IN) :: fsource
        CHARACTER(LEN=len_filename), INTENT(IN) :: ftarget
        CHARACTER(LEN=len_syscommand) :: mycommand, CMSG
        INTEGER :: EXSTAT, CSTAT, ICHK

        ICHK = INDEX(TRIM(fsource),TRIM(HEADDIR))
        IF(ICHK .EQ. 0) THEN
            WRITE(*,*) "ERROR) FMOVE REQUESTED FOR PATH OUTSIDE HEADDIR"
            WRITE(*,*) "ERROR) HEADDIR : ", TRIM(HEADDIR)
            STOP
        END IF

        ICHK = INDEX(TRIM(ftarget),TRIM(HEADDIR))
        IF(ICHK .EQ. 0) THEN
            WRITE(*,*) "ERROR) FMOVE REQUESTED FOR PATH OUTSIDE HEADDIR"
            WRITE(*,*) "ERROR) HEADDIR : ", TRIM(HEADDIR)
            STOP
        END IF

        mycommand = 'mv '//TRIM(ADJUSTL(fsource))//' '//TRIM(ADJUSTL(ftarget))

        CALL EXECUTE_COMMAND_LINE(TRIM(mycommand),EXITSTAT=EXSTAT,CMDSTAT=CSTAT,CMDMSG=CMSG)

        IF(CSTAT .NE. 0) THEN
            WRITE(*,*) "ERROR) COULD NOT EXECUTE FMOVE COMMAND.."
            WRITE(*,*) "ERROR) COMMAND  : ", TRIM(mycommand)
            WRITE(*,*) "ERROR) EXIT MSG : ", TRIM(CMSG)
            STOP
        END IF
        IF(EXSTAT .NE. 0) THEN
            WRITE(*,*) "ERROR) EXECUTED MOVE COMMAND THROWS AN ERROR.."
            WRITE(*,*) "ERROR) COMMAND  : ", TRIM(mycommand)
            STOP
        END IF
    END SUBROUTINE FMOVE

    SUBROUTINE DMOVE(dsource,dtarget)
        CHARACTER(LEN=len_filename), INTENT(IN) :: dsource
        CHARACTER(LEN=len_filename), INTENT(IN) :: dtarget
        CHARACTER(LEN=len_syscommand) :: mycommand, CMSG
        INTEGER :: EXSTAT, CSTAT, ICHK

        ICHK = INDEX(TRIM(dsource),TRIM(HEADDIR))
        IF(ICHK .EQ. 0) THEN
            WRITE(*,*) "ERROR) DMOVE REQUESTED FOR PATH OUTSIDE HEADDIR"
            WRITE(*,*) "ERROR) HEADDIR : ", TRIM(HEADDIR)
            STOP
        END IF

        ICHK = INDEX(TRIM(dtarget),TRIM(HEADDIR))
        IF(ICHK .EQ. 0) THEN
            WRITE(*,*) "ERROR) DMOVE REQUESTED FOR PATH OUTSIDE HEADDIR"
            WRITE(*,*) "ERROR) HEADDIR : ", TRIM(HEADDIR)
            STOP
        END IF

        mycommand = 'mv '//TRIM(ADJUSTL(dsource))//' '//TRIM(ADJUSTL(dtarget))

        CALL EXECUTE_COMMAND_LINE(TRIM(mycommand),EXITSTAT=EXSTAT,CMDSTAT=CSTAT,CMDMSG=CMSG)

        IF(CSTAT .NE. 0) THEN
            WRITE(*,*) "ERROR) COULD NOT EXECUTE DMOVE COMMAND.."
            WRITE(*,*) "ERROR) COMMAND  : ", TRIM(mycommand)
            WRITE(*,*) "ERROR) EXIT MSG : ", TRIM(CMSG)
            STOP
        END IF

        IF(EXSTAT .NE. 0) THEN
            WRITE(*,*) "ERROR) EXECUTED DMOVE COMMAND THROWS AN ERROR.."
            WRITE(*,*) "ERROR) COMMAND  : ", TRIM(mycommand)
            STOP
        END IF
    END SUBROUTINE DMOVE

    SUBROUTINE CREATEDIR(dir)
        CHARACTER(LEN=len_filename), INTENT(IN) :: dir
        CHARACTER(LEN=len_syscommand) :: mycommand, CMSG
        INTEGER :: EXSTAT, CSTAT, ICHK

        ICHK = INDEX(TRIM(dir),TRIM(HEADDIR))
        IF(ICHK .EQ. 0) THEN
            WRITE(*,*) "ERROR) CREATEDIR REQUESTED FOR PATH OUTSIDE HEADDIR"
            WRITE(*,*) "ERROR) HEADDIR : ", TRIM(HEADDIR)
            STOP
        END IF

        mycommand = 'mkdir -p '//TRIM(ADJUSTL(dir))

        CALL EXECUTE_COMMAND_LINE(TRIM(mycommand),EXITSTAT=EXSTAT,CMDSTAT=CSTAT,CMDMSG=CMSG)

        IF(CSTAT .NE. 0) THEN
            WRITE(*,*) "ERROR) COULD NOT EXECUTE CREATEDIR COMMAND.."
            WRITE(*,*) "ERROR) COMMAND  : ", TRIM(mycommand)
            WRITE(*,*) "ERROR) EXIT MSG : ", TRIM(CMSG)
            STOP
        END IF

        IF(EXSTAT .NE. 0) THEN
            WRITE(*,*) "ERROR) EXECUTED MKDIR COMMAND THROWS AN ERROR.."
            WRITE(*,*) "ERROR) COMMAND  : ", TRIM(mycommand)
            STOP
        END IF
    END SUBROUTINE CREATEDIR
END MODULE sysutils