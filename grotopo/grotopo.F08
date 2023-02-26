
! A MODULE TO READ AND WRITE TOPOLOGY INFORMATION (IN GROMACS FORMAT) RELATED TO ATOM TYPES
! IT ALSO POPULATES THE LJPAR DATA STRUCTURE

! STARTED ON 24 JUNE 2019
! UPDATED ON 05 JULY 2019

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! THE API !!!!!
    ! TYPE(grotype), ALLOCATABLE :: itppars(:)
    ! ALLOCATE(itppars(numpars))
    
    ! CALL readgtp('/home/nikhil/projects/2_diether/codes/tunelj/datafiles/guesspar.itp',itppars)
    ! CALL writegtp('/home/nikhil/projects/2_diether/codes/tunelj/grotopo/testout.itp',itppars)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !**** GTP FILE FORMAT ****! A CUSTOM DEFINED FORMAT
    !; name  at.num  mass  charge  ptype  sigma  epsilon    lsigma usigma ssigma    lepsilon uepsilon sepsilon

    ! THE TEMPLATE TOPOLOGY FILE USES INCLUDE STATMENT -- NOW SPECIFIED IN GLOBAL MODULE
    ! #include "par.gtp"

MODULE grotopo
USE myglobal
IMPLICIT NONE

    TYPE :: grotype
        TYPE(ljpartype)            :: gropar
        CHARACTER(LEN=len_partype) :: gparticle
        REAL(KIND=kind_real)       :: gmass
        REAL(KIND=kind_real)       :: gcharge
        INTEGER                    :: gelnum
    END TYPE grotype

    CONTAINS
    ! READS FROM THE INPUT GTP FILE AND LOADS THE DATA INTO DATA STRUCTURE
    SUBROUTINE READGTP(filename,itppars)
        CHARACTER(LEN=len_filename), INTENT(IN) ::filename
        TYPE(grotype),INTENT(INOUT) :: itppars(:)
        INTEGER :: ITPUNIT, ios, i

        OPEN(NEWUNIT=ITPUNIT,FILE=TRIM(filename),IOSTAT=ios,STATUS='OLD',ACTION='READ')
        IF(ios .NE. 0) THEN
            WRITE(*,'(A,A)') "ERROR READING FROM ", TRIM(filename)
            STOP
        END IF

        DO i = 1, SIZE(itppars)
            READ(ITPUNIT,*,IOSTAT=ios) itppars(i)%gropar%parname, &
            itppars(i)%gelnum, itppars(i)%gmass, itppars(i)%gcharge, itppars(i)%gparticle, &
            itppars(i)%gropar%par(1), itppars(i)%gropar%par(2), &
            itppars(i)%gropar%lwrbnd(1), itppars(i)%gropar%uprbnd(1), itppars(i)%gropar%sens(1), & 
            itppars(i)%gropar%lwrbnd(2), itppars(i)%gropar%uprbnd(2), itppars(i)%gropar%sens(2)
            IF(ios .NE. 0) THEN
                IF(IS_IOSTAT_END(ios)) THEN
                    WRITE(*,'(A,A)') "ERROR) UNABLE TO READ FROM FILE :", TRIM(filename)
                    WRITE(*,'(A,I0,A,I0)') "ERROR) END OF FILE REACHED AT:", i," BEFORE READING: ", &
                                           SIZE(itppars)," SPECIFIED PARAMETERS" 
                    STOP
                END IF
                WRITE(*,'(A,A)') "ERROR) UNABLE TO READ FROM FILE :", TRIM(filename)
                STOP
            END IF
        END DO
        CLOSE(ITPUNIT)
    END SUBROUTINE READGTP

    ! WRITES THE LJ PARAMETERS TO A TEMPORARY GTP FILE -- INCLUDES THE BOUNDS AND THE SENSITIVITY
    SUBROUTINE WRITEGTP(filename,itppars)
        CHARACTER(LEN=len_filename), INTENT(IN) ::filename
        TYPE(grotype),INTENT(INOUT) :: itppars(:)
        INTEGER :: GTPUNIT, ios, i

        OPEN(NEWUNIT=GTPUNIT,FILE=TRIM(filename),IOSTAT=ios)
        IF(ios .NE. 0) THEN
            WRITE(*,*) "ERROR OPENING ", TRIM(filename)
            STOP
        END IF
        119     FORMAT(A4,2X,I4,2X,2(F8.4,2X),A1,2X,F8.4,2X,F10.6,6X,6(F10.6,2X))
        DO i = 1, SIZE(itppars)
            WRITE(GTPUNIT,119) TRIM(ADJUSTL(itppars(i)%gropar%parname)), &
                               itppars(i)%gelnum, itppars(i)%gmass, itppars(i)%gcharge, &
                               TRIM(ADJUSTL(itppars(i)%gparticle)), &
                               itppars(i)%gropar%par(1), itppars(i)%gropar%par(2), &
                               itppars(i)%gropar%lwrbnd(1), itppars(i)%gropar%uprbnd(1), itppars(i)%gropar%sens(1), & 
                               itppars(i)%gropar%lwrbnd(2), itppars(i)%gropar%uprbnd(2), itppars(i)%gropar%sens(2)
        END DO
        CLOSE(GTPUNIT)
    END SUBROUTINE WRITEGTP

    ! WRITES THE LJ PARAMETERS TO A TEMPORARY ITP FILE TO BE USED TO CREATE FULL TOPOLOGY FILE
    SUBROUTINE WRITEITP(filename,itppars)
        CHARACTER(LEN=len_filename), INTENT(IN) ::filename
        TYPE(grotype),INTENT(INOUT) :: itppars(:)
        INTEGER :: ITPUNIT, ios, i

        OPEN(NEWUNIT=ITPUNIT,FILE=TRIM(filename),IOSTAT=ios)
        IF(ios .NE. 0) THEN
            WRITE(*,*) "ERROR OPENING ", TRIM(filename)
            STOP
        END IF
        109     FORMAT(A4,2X,I4,2X,2(F8.4,2X),A1,2X,F8.4,2X,F10.6)
        DO i = 1, SIZE(itppars)
            WRITE(ITPUNIT,109) TRIM(ADJUSTL(itppars(i)%gropar%parname)), &
                               itppars(i)%gelnum, itppars(i)%gmass, itppars(i)%gcharge, &
                               TRIM(ADJUSTL(itppars(i)%gparticle)), &
                               itppars(i)%gropar%par(1), itppars(i)%gropar%par(2)
        END DO
        CLOSE(ITPUNIT)
    END SUBROUTINE WRITEITP

    SUBROUTINE READENERGY(filename, myenergy)
        CHARACTER(LEN=len_filename), INTENT(IN) ::filename
        REAL(KIND=kind_real), INTENT(INOUT)     :: myenergy(:)
        
        INTEGER :: ENERUNIT, ios, i
        CHARACTER(LEN=len_record) :: thisline, tempstr

        OPEN(NEWUNIT=ENERUNIT,FILE=TRIM(filename),IOSTAT=ios)
        IF(ios .NE. 0) THEN
            WRITE(*,*) "ERROR OPENING ", TRIM(filename)
            STOP
        END IF

        i = 0
        DO 
            READ(ENERUNIT,'(A)',IOSTAT=ios) thisline
            
            ! CYCLE BLANK LINES
            IF(LEN_TRIM(thisline) .EQ. 0 ) CYCLE

            ! CYCLE COMMENT LINES AND CHECK FOR REAL NUMBERS
            tempstr = ADJUSTL(thisline)
            SELECT CASE(tempstr(1:1))
                CASE('#','@','!')
                    ! DEBUGL i = i + 1
                    CYCLE
                CASE('0','1','2','3','4','5','6','7','8','9','.')
                    ! DEBUGL i = i + 1
                    EXIT
                CASE DEFAULT
                    WRITE(*,*) "ERROR) UNKNOWN STRING IN .XVG FILE "//tempstr(1:1)
                    WRITE(*,*) "ERROR) FILE :"//TRIM(filename)
                    STOP
            END SELECT

            ! IF FILE ENDS BEFORE WE GET SOME DATA THROW AND ERROR AND STOP
            IF(IS_IOSTAT_END(ios)) THEN
                WRITE(*,*) "ERROR) COULDNOT FIND A VALID ENERGY LINE IN .XVG FILE"
                WRITE(*,*) "ERROR) FILE :"//TRIM(filename)
                STOP
            END IF
        END DO
        ! DEBUGL WRITE(*,'(A,I0)') "NUMBER OF LINES IN XVG FILE.. ", i
        ! READ THE VALUES INTO AN ARRAY
        READ(thisline,*,IOSTAT=ios) myenergy
        IF(ios .NE. 0) THEN
            WRITE(*,*) "ERROR) .XVG ENERGY FILE IS POSSIBLY CORRUPT.."
            WRITE(*,*) "ERROR) FILE :"//TRIM(filename)
            STOP
        END IF
        CLOSE(ENERUNIT)
    END SUBROUTINE READENERGY

    SUBROUTINE CHECKPAR(pars,LOGUNIT)
        TYPE(ljpartype),INTENT(IN) :: pars(:)
        INTEGER, INTENT(IN) :: LOGUNIT
        INTEGER :: i, numlj

        numlj = SIZE(pars)

        DO i = 1, numlj
            !ERRORS
            IF((pars(i)%uprbnd(1) .LT. pars(i)%lwrbnd(1)) .OR. (pars(i)%uprbnd(2) .LT. pars(i)%lwrbnd(2))) THEN
                WRITE(*,'(A,A)') " CHECKPAR) ERROR : LOWER BOUND IS GREATER THAN UPPER BOUND. PAR: ", TRIM(pars(i)%parname) 
                WRITE(LOGUNIT,'(A,A)') " CHECKPAR) ERROR : LOWER BOUND IS GREATER THAN UPPER BOUND. PAR: ", TRIM(pars(i)%parname) 
                STOP
            ELSE IF((pars(i)%par(1) .LT. pars(i)%lwrbnd(1)) .OR. (pars(i)%par(1) .GT. pars(i)%uprbnd(1))) THEN
                WRITE(*,'(A,A)') " CHECKPAR) ERROR : PARAMETER OUT OF BOUNDS. PAR: ", TRIM(pars(i)%parname) 
                WRITE(LOGUNIT,'(A,A)') " CHECKPAR) ERROR : PARAMETER OUT OF BOUNDS. PAR: ", TRIM(pars(i)%parname) 
                STOP
            ELSE IF((pars(i)%par(2) .LT. pars(i)%lwrbnd(2)) .OR. (pars(i)%par(2) .GT. pars(i)%uprbnd(2))) THEN
                WRITE(*,'(A,A)') " CHECKPAR) ERROR : PARAMETER OUT OF BOUNDS. PAR: ", TRIM(pars(i)%parname) 
                WRITE(LOGUNIT,'(A,A)') " CHECKPAR) ERROR : PARAMETER OUT OF BOUNDS. PAR: ", TRIM(pars(i)%parname) 
                STOP
            ELSE IF(pars(i)%sens(1) .GT. (pars(i)%uprbnd(1)-pars(i)%lwrbnd(1))) THEN
                WRITE(*,'(A,A)') " CHECKPAR) ERROR : SENSITIVITY GREATER THAN RANGE. PAR: ", TRIM(pars(i)%parname) 
                WRITE(LOGUNIT,'(A,A)') " CHECKPAR) ERROR : SENSITIVITY GREATER THAN RANGE. PAR: ", TRIM(pars(i)%parname) 
                STOP
            ELSE IF(pars(i)%sens(2) .GT. (pars(i)%uprbnd(2)-pars(i)%lwrbnd(2))) THEN
                WRITE(*,'(A,A)') " CHECKPAR) ERROR : SENSITIVITY GREATER THAN RANGE. PAR: ", TRIM(pars(i)%parname) 
                WRITE(LOGUNIT,'(A,A)') " CHECKPAR) ERROR : SENSITIVITY GREATER THAN RANGE. PAR: ", TRIM(pars(i)%parname) 
                STOP
            END IF
            ! WARNINGS
            IF(pars(i)%sens(1) .GT. 0.5D0*(pars(i)%uprbnd(1)-pars(i)%lwrbnd(1))) THEN
                WRITE(*,'(A,A)') " CHECKPAR) WARNING!! : SENSITIVITY GREATER THAN HALF RANGE. PAR: ", TRIM(pars(i)%parname) 
                WRITE(LOGUNIT,'(A,A)') " CHECKPAR) WARNING!! : SENSITIVITY GREATER THAN HALF RANGE. PAR: ", TRIM(pars(i)%parname) 
            ENDIF
            IF(pars(i)%sens(2) .GT. 0.5D0*(pars(i)%uprbnd(2)-pars(i)%lwrbnd(2))) THEN
                WRITE(*,'(A,A)') " CHECKPAR) WARNING!! : SENSITIVITY GREATER THAN HALF RANGE. PAR: ", TRIM(pars(i)%parname) 
                WRITE(LOGUNIT,'(A,A)') " CHECKPAR) WARNING!! : SENSITIVITY GREATER THAN HALF RANGE. PAR: ", TRIM(pars(i)%parname) 
            ENDIF
        END DO
    END SUBROUTINE CHECKPAR

    SUBROUTINE GETPARARR(myljpar, mypararr)
        TYPE(ljpartype), INTENT(IN)       :: myljpar(:)
        REAL(KIND=kind_real), INTENT(OUT) :: mypararr(:)
        INTEGER :: i, ljparsize, pararrsize

        ljparsize = SIZE(myljpar) ; pararrsize = SIZE(mypararr)
        IF(MOD(pararrsize,2) .NE. 0) THEN
            WRITE(*,'(A)')"GETPARRARR) ERROR : SIZE OF THE PAR ARRAY RECEIVED IS NOT A MULTIPLE OF TWO"
            STOP
        END IF
        IF(pararrsize/2 .NE. ljparsize) THEN
            WRITE(*,'(A)')"GETPARRARR) ERROR : SIZE OF THE PAR ARRAY AND LJPAR ARE NOT COMPATIBLE"
            STOP
        END IF

        DO i = 1, ljparsize
            mypararr(2*i-1) = myljpar(i)%par(1)
            mypararr(2*i)   = myljpar(i)%par(2)
        END DO
    END SUBROUTINE GETPARARR

    SUBROUTINE PUTPARARR(mypararr, myljpar)
        REAL(KIND=kind_real), INTENT(IN) :: mypararr(:)
        TYPE(ljpartype), INTENT(OUT)     :: myljpar(:)
        INTEGER :: i, pararrsize

        pararrsize = SIZE(mypararr)
        IF(MOD(pararrsize,2) .NE. 0) THEN
            WRITE(*,'(A)')"PUTPARRARR) ERROR : SIZE OF THE PAR ARRAY RECEIVED IS NOT A MULTIPLE OF TWO"
            STOP
        END IF
        IF(pararrsize/2 .NE. SIZE(myljpar)) THEN
            WRITE(*,'(A)')"PUTPARRARR) ERROR : SIZE OF THE PAR ARRAY AND LJPAR ARE NOT COMPATIBLE"
            STOP
        END IF

        DO i = 1, pararrsize/2
            myljpar(i)%par(1) = mypararr(2*i-1)
            myljpar(i)%par(2) = mypararr(2*i)
        END DO
    END SUBROUTINE PUTPARARR

    SUBROUTINE GETSENSARR(myljpar, mysensarr)
        TYPE(ljpartype), INTENT(IN)       :: myljpar(:)
        REAL(KIND=kind_real), INTENT(OUT) :: mysensarr(:)
        INTEGER :: i, ljparsize, sensarrsize

        ljparsize = SIZE(myljpar) ; sensarrsize = SIZE(mysensarr)
        IF(MOD(sensarrsize,2) .NE. 0) THEN
            WRITE(*,'(A)')"GETSENSRARR) ERROR : SIZE OF THE SENS ARRAY RECEIVED IS NOT A MULTIPLE OF TWO"
            STOP
        END IF
        IF(sensarrsize/2 .NE. ljparsize) THEN
            WRITE(*,'(A)')"GETSENSRARR) ERROR : SIZE OF THE SENS ARRAY AND LJPAR ARE NOT COMPATIBLE"
            STOP
        END IF

        DO i = 1, ljparsize
            mysensarr(2*i-1) = myljpar(i)%sens(1)
            mysensarr(2*i)   = myljpar(i)%sens(2)
        END DO
    END SUBROUTINE GETSENSARR

    SUBROUTINE PUTSENSARR(mysensarr, myljpar)
        REAL(KIND=kind_real), INTENT(IN) :: mysensarr(:)
        TYPE(ljpartype), INTENT(OUT)     :: myljpar(:)
        INTEGER :: i, sensarrsize

        sensarrsize = SIZE(mysensarr)
        IF(MOD(sensarrsize,2) .NE. 0) THEN
            WRITE(*,'(A)')"PUTSENSARR) ERROR : SIZE OF THE SENS ARRAY RECEIVED IS NOT A MULTIPLE OF TWO"
            STOP
        END IF
        IF(sensarrsize/2 .NE. SIZE(myljpar)) THEN
            WRITE(*,'(A)')"PUTSENSARR) ERROR : SIZE OF THE SENS ARRAY AND LJPAR ARE NOT COMPATIBLE"
            STOP
        END IF

        DO i = 1, sensarrsize/2
            myljpar(i)%sens(1) = mysensarr(2*i-1)
            myljpar(i)%sens(2) = mysensarr(2*i)
        END DO
    END SUBROUTINE PUTSENSARR
END MODULE grotopo