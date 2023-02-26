
! PROGRAM TO TEST GROTOPO MODULE

PROGRAM test
USE myglobal
USE grotopo
IMPLICIT NONE

    INTEGER, PARAMETER :: numpars = 2

    TYPE(grotype), ALLOCATABLE :: itppars(:)

    ALLOCATE(itppars(numpars))

    CALL readitp('/home/nikhil/projects/2_diether/codes/tunelj/datafiles/guesspar.itp',itppars)
    CALL writeitp('/home/nikhil/projects/2_diether/codes/tunelj/grotopo/testout.itp',itppars)

END PROGRAM test