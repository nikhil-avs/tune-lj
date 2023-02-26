
! GLOBAL MODULE TO CONTAIN GLOBAL VARIABLES, CONSTANTS AND PARAMETERS

! STARTED ON 24 JUNE 2019
! UPDATED ON 05 JULY 2019

MODULE myglobal
USE ISO_FORTRAN_ENV
IMPLICIT NONE
    INTEGER, PARAMETER          :: len_partype    = 1
    INTEGER, PARAMETER          :: len_parname    = 10
    INTEGER, PARAMETER          :: len_scanname   = 20
    INTEGER, PARAMETER          :: len_format     = 100
    INTEGER, PARAMETER          :: len_dirname    = 400
    INTEGER, PARAMETER          :: len_filename   = 400
    INTEGER, PARAMETER          :: len_record     = 400
    INTEGER, PARAMETER          :: len_syscommand = 4000
    INTEGER, PARAMETER          :: kind_real      = 8
    INTEGER, PARAMETER          :: maxflags       = 2
    INTEGER, PARAMETER          :: ENERGYPARTS    = 9
    INTEGER, PARAMETER          :: IBONDENERGY    = 2
    INTEGER, PARAMETER          :: IANGLEENERGY   = 3
    INTEGER, PARAMETER          :: IDIHEDENERGY   = 4
    ! INTEGER, PARAMETER          :: ILJENERGY      = 4
    ! INTEGER, PARAMETER          :: ICOULENERGY    = 5
    INTEGER, PARAMETER          :: IPOTENERGY     = 9
    CHARACTER(LEN=4), PARAMETER :: CONF_EXT       = '.gro'
    CHARACTER(LEN=7), PARAMETER :: INCLUDE_GTP_FILE = 'par.gtp'
    INTEGER, PARAMETER          :: max_move_tries   = 2000

    TYPE :: ljpartype
        CHARACTER(LEN=len_parname) :: parname
        CHARACTER(LEN=len_partype) :: partype
        REAL(KIND=kind_real)       :: par(2)
        REAL(KIND=kind_real)       :: lwrbnd(2)
        REAL(KIND=kind_real)       :: uprbnd(2)
        REAL(KIND=kind_real)       :: sens(2)
    END TYPE ljpartype

    TYPE :: conftype
        CHARACTER(LEN=len_filename)              :: cscanname
        CHARACTER(LEN=len_dirname)               :: cdirname
        CHARACTER(LEN=len_filename)              :: cfilename
        INTEGER                                  :: numconfs
        CHARACTER(LEN=len_filename), ALLOCATABLE :: confname(:)
        REAL(KIND=kind_real), ALLOCATABLE        :: cdist(:)
        CHARACTER(LEN=len_filename), ALLOCATABLE :: topolname(:)
        CHARACTER(LEN=len_filename), ALLOCATABLE :: conffilename(:)
    END TYPE conftype

    TYPE :: qmdatatype
        CHARACTER(LEN=len_scanname)       :: qscanname
        CHARACTER(LEN=len_filename)       :: qfilename
        INTEGER                           :: numdata
        REAL(KIND=kind_real), ALLOCATABLE :: qdist(:)
        REAL(KIND=kind_real), ALLOCATABLE :: qenergy(:)
        REAL(KIND=kind_real), ALLOCATABLE :: qweight(:)
    END TYPE qmdatatype

    TYPE :: betatype
        INTEGER                           :: numpoints
        INTEGER                           :: currpnt
        INTEGER, ALLOCATABLE              :: biter(:)
        REAL(KIND=kind_real), ALLOCATABLE :: beta(:)
        REAL(KIND=kind_real)              :: binit
        REAL(KIND=kind_real)              :: bcurr
        REAL(KIND=kind_real)              :: brate
    END TYPE betatype

    CHARACTER(len=len_dirname)  :: HEADDIR, RUNDIR
    CHARACTER(len=len_filename) :: ALLSPACE_FILE
    CHARACTER(len=len_dirname)  :: ALLSPACE_DIR
    CHARACTER(len=len_scanname) :: ALLSPACE_SCAN 
    CHARACTER(len=len_record)   :: ALLSPACE_REC 
    CHARACTER(len=len_dirname)  :: GROMACSEXE
    CHARACTER(LEN=1)            :: mynewline
END MODULE myglobal

