
! MAIN PROGRAM WHERE THE STUFF HAPPENS!

! STARTED ON 24 JUNE 2019
! UPDATED ON 05 JULY 2019

PROGRAM main
USE myglobal
USE sysutils
USE my_cl_parser
USE grotopo
IMPLICIT NONE
    ! FILE AND FOLDER NAMES
    CHARACTER(len=len_filename), DIMENSION(maxflags) :: strinp
    CHARACTER(len=len_filename)  :: MDHEADDIR, QMHEADDIR, SENSDIR, DATADIR
    CHARACTER(len=len_filename)  :: qm_info_file, conf_info_file
    CHARACTER(len=len_filename)  :: mdpfilename, gtpfilename, currgtpfilename
    CHARACTER(len=len_filename)  :: tempfile
    INTEGER                      :: LOGUNIT

    ! PARAMETER RELATED DATA STRUCTURES
    INTEGER                      :: NUMLJPAR
    TYPE(grotype), ALLOCATABLE   :: firstitppar(:), curritppar(:), nextitppar(:)
    TYPE(ljpartype), ALLOCATABLE :: firstpar(:), currpar(:), nextpar(:), minpar(:)
    REAL(KIND=kind_real), ALLOCATABLE :: pararr(:), sensarr(:)
    ! This is to track the minimum RMSE paramter -- miniter(1) = 0 (equil) = 1 (prod)
    INTEGER                           :: miniter(2)

    ! QM AND MD ENERGY DATA STRUCTURES
    INTEGER                      :: NUMSCANS
    TYPE(qmdatatype),ALLOCATABLE :: qmdata(:), mdfirst(:), mdnow(:), mdnext(:)

    ! CONFIGURATIONS RELATED DATA STRUCTURES
    TYPE(conftype),ALLOCATABLE   :: myconfs(:)

    ! MONE CARLO RELATED VARIABLES
    REAL(KIND=kind_real) :: firstcumerr, currcumerr, nextcumerr, delcumerr, mincumerr
    TYPE(betatype)       :: mcbeta
    INTEGER              :: mciter, MCEQUIL, MCPROD
    INTEGER              :: mcaccept, mcreject, mc_anneal_num
    LOGICAL              :: ACCEPTED, MIN_PROD_START
    INTEGER              :: IACC

    ! ACCEPTANCE RATIO VARIABLES
    REAL(KIND=kind_real) :: curr_acc_ratio, prev_acc_ratio, ref_acc_ratio
    REAL(KIND=kind_real) :: err_acc_ratio, interr_acc_ratio, diferr_acc_ratio, preverr_acc_ratio
    REAL(KIND=kind_real) :: kp, ki, kd, corr_factor, sens_scale
    INTEGER              :: curr_acc_iter
    INTEGER              :: CONTROL_ACC_EVERY, UNIT_ACCLOG, UNIT_SENSLOG
    REAL(KIND=kind_real) :: MINSENS_SIG = 0.00001D0, MINSENS_EPS = 0.00001D0
    LOGICAL              :: CONTROL_ACC

    REAL(KIND=kind_real) :: myrandnum(2)
    INTEGER :: MCALLENER, MCALLPAR, MCACCENER, MCACCPAR, MCREJENER, MCREJPAR
    INTEGER :: IMC_EQUIL=0, IMC_PROD=1

    ! BOOK KEEPING VARIABLES
    INTEGER :: i, j, ios, mysecs
    INTEGER, ALLOCATABLE :: myrandseed(:)
    CHARACTER(LEN=len_format) :: myformat, myformat2

    ! GLOBAL INITIALIZATION
    ALLSPACE_FILE = REPEAT(' ',len_filename)
    ALLSPACE_SCAN = REPEAT(' ',len_scanname)
    ALLSPACE_DIR  = REPEAT(' ',len_dirname)
    ALLSPACE_REC  = REPEAT(' ',len_record)
    mynewline = NEW_LINE(HEADDIR)

    CALL SYSTEM_CLOCK(COUNT=mysecs)
    CALL RANDOM_SEED(SIZE=i)
    ALLOCATE(myrandseed(i)) ; myrandseed = 0 ; myrandseed = myrandseed + mysecs
    CALL RANDOM_SEED(PUT=myrandseed)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!************************************************* INITIALIZATION PART *************************************************!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    WRITE(*,'(A)') "INIT) INITIALIZATION STARTS.."
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! READ AND STORE THE COMMAND LINE FLAGS
    CALL INPUT_FLAGS(strinp)
    OPEN(NEWUNIT=LOGUNIT,FILE=TRIM(strinp(2))//".log",IOSTAT=ios)
    IF(ios .NE. 0) THEN
        WRITE(*,'(A,A)')       "ERROR) CANNOT OPEN LOG FILE : ", TRIM(strinp(2))
        WRITE(LOGUNIT,'(A,A)') "ERROR) CANNOT OPEN LOG FILE : ", TRIM(strinp(2))
        STOP
    END IF
    
    WRITE(*,'(A,A)')       "INIT) READING JOB CONTROL VALUES FROM : ", TRIM(strinp(1))
    WRITE(LOGUNIT,'(A,A)') "INIT) READING JOB CONTROL VALUES FROM : ", TRIM(strinp(1))
    
    CALL PARSEINP(strinp(1))
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!********************************** SETTING UP THE WORKSPACE *******************************!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    RUNDIR = TRIM(HEADDIR)//'/'//TRIM(strinp(2))//'/'
    SENSDIR = TRIM(RUNDIR)//'sens/'
    DATADIR = TRIM(RUNDIR)//'data/'
    CALL CREATEDIR(RUNDIR)
    WRITE(*,'(A)')       "INIT) CREATED WORKSPACE.. "
    WRITE(LOGUNIT,'(A)') "INIT) CREATED WORKSPACE.. "
    CALL CREATEDIR(SENSDIR)
    WRITE(*,'(A)')       "INIT) CREATED SENSDIR.. "
    WRITE(LOGUNIT,'(A)') "INIT) CREATED SENSDIR.. "
    CALL CREATEDIR(DATADIR)
    WRITE(*,'(A)')       "INIT) CREATED DATADIR.. "
    WRITE(LOGUNIT,'(A)') "INIT) CREATED DATADIR.. "
    !*********************************************************************************************!
    CALL PRINT_LOG_PARAMS()

    CALL SORTBETA()
    WRITE(*,'(A)')       "MCPROD) SORTED AND INITIALIZED BETA ARRAY"
    WRITE(LOGUNIT,'(A)') "MCPROD) SORTED AND INITIALIZED BETA ARRAY"

    WRITE(*,'(A)')       "--> ITERATION    BETA"
    WRITE(LOGUNIT,'(A)') "--> ITERATION    BETA"
    DO i = 0, mcbeta%numpoints+1
        WRITE(*,'(A,I4,2X,F14.6)')       "--> ", mcbeta%biter(i), mcbeta%beta(i)
        WRITE(LOGUNIT,'(A,I4,2X,F14.6)') "--> ", mcbeta%biter(i), mcbeta%beta(i)
    END DO

    !*********************************************************************************************!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!************************************* PARAMETER RELATED ************************************!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! INITIALIZE THE ITP PARAMATERS DATA STRUCTURE AND LOAD THE RELEVANT INFORMATION INTO PARAMETER VECTOR
    WRITE(*,'(A)') "INIT) READING GUESS PARAMETERS FROM GTP FILE : ", TRIM(MDHEADDIR)//'/'//TRIM(gtpfilename)
    WRITE(LOGUNIT,'(A)') "INIT) READING GUESS PARAMETERS FROM GTP FILE : ", TRIM(MDHEADDIR)//'/'//TRIM(gtpfilename)
    
    ALLOCATE(firstitppar(NUMLJPAR), curritppar(NUMLJPAR), nextitppar(NUMLJPAR))
    ALLOCATE(firstpar(NUMLJPAR), currpar(NUMLJPAR), nextpar(NUMLJPAR), minpar(NUMLJPAR))
    ALLOCATE(pararr(2*NUMLJPAR),sensarr(2*NUMLJPAR))
    
    tempfile = ALLSPACE_FILE
    tempfile = TRIM(MDHEADDIR)//'/'//TRIM(gtpfilename)
    CALL READGTP(tempfile,firstitppar)
    
    curritppar(:) = firstitppar(:)
    nextitppar(:) = firstitppar(:)
    DO i = 1, NUMLJPAR
        firstpar(i) = firstitppar(i)%gropar
        currpar(i) = firstitppar(i)%gropar
        nextpar(i) = firstitppar(i)%gropar
    END DO

    CALL CHECKPAR(firstpar,LOGUNIT)
    WRITE(*,'(A)')       "INIT) CONSISTANCY OF PARAMETER BOUNDS CHECKED.. ALL OK!"
    WRITE(LOGUNIT,'(A)') "INIT) CONSISTANCY OF PARAMETER BOUNDS CHECKED.. ALL OK!"

    WRITE(*,'(A,A)')       "INIT) !!!!!!!!!!!!! INITIAL GUESS PARAMTERS !!!!!!!!!!!!!!!!!! "
    WRITE(LOGUNIT,'(A,A)') "INIT) !!!!!!!!!!!!! INITIAL GUESS PARAMTERS !!!!!!!!!!!!!!!!!! "

    WRITE(*,'(A)')       "--> PAR     SIGMA         LOWER        UPPER         SENS"
    WRITE(LOGUNIT,'(A)') "--> PAR     SIGMA         LOWER        UPPER         SENS"
    DO i = 1, NUMLJPAR
        WRITE(*,'(A,A,4(F12.6,2X))')       "--> ", TRIM(firstitppar(i)%gropar%parname),&
                        firstitppar(i)%gropar%par(1), firstitppar(i)%gropar%lwrbnd(1),&
                        firstitppar(i)%gropar%uprbnd(1), firstitppar(i)%gropar%sens(1)
        WRITE(LOGUNIT,'(A,A,4(F12.6,2X)))') "--> ", TRIM(firstitppar(i)%gropar%parname),&
                        firstitppar(i)%gropar%par(1), firstitppar(i)%gropar%lwrbnd(1),&
                        firstitppar(i)%gropar%uprbnd(1), firstitppar(i)%gropar%sens(1)
    END DO

    WRITE(*,'(A)')       "--> PAR    EPSILON        LOWER        UPPER         SENS"
    WRITE(LOGUNIT,'(A)') "--> PAR    EPSILON        LOWER        UPPER         SENS"
    DO i = 1, NUMLJPAR
        WRITE(*,'(A,A,4(F12.6,2X))')       "--> ", TRIM(firstitppar(i)%gropar%parname),&
                        firstitppar(i)%gropar%par(2), firstitppar(i)%gropar%lwrbnd(2),&
                        firstitppar(i)%gropar%uprbnd(2), firstitppar(i)%gropar%sens(2)
        WRITE(LOGUNIT,'(A,A,4(F12.6,2X)))') "--> ", TRIM(firstitppar(i)%gropar%parname),&
                        firstitppar(i)%gropar%par(2), firstitppar(i)%gropar%lwrbnd(2),&
                        firstitppar(i)%gropar%uprbnd(2), firstitppar(i)%gropar%sens(2)
    END DO
    WRITE(*,'(A,A)')       "INIT) ************* INITIAL GUESS PARAMTERS ****************** "
    WRITE(LOGUNIT,'(A,A)') "INIT) ************* INITIAL GUESS PARAMTERS ****************** "

    !*********************************************************************************************!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!************************************* QM DATA RELATED *************************************!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! THE ORDER OF CONFS AND QM CAN BE DIFFERENT AND CAN BE UNSORTED
    tempfile = ALLSPACE_FILE
    tempfile = TRIM(QMHEADDIR)//'/'//TRIM(qm_info_file)

    WRITE(*,'(A,A)') "INIT) READING QM INFO FILE : ", TRIM(tempfile)
    WRITE(LOGUNIT,'(A,A)') "INIT) READING QM INFO FILE : ", TRIM(tempfile)

    ALLOCATE(qmdata(NUMSCANS),mdfirst(NUMSCANS))
    ALLOCATE(mdnow(NUMSCANS),mdnext(NUMSCANS))
    ALLOCATE(myconfs(NUMSCANS))
    
    CALL READQMDATA(tempfile,qmdata)
    CALL SORTQMDATA(qmdata)

    DO i = 1, NUMSCANS
        WRITE(*,'(A,I0,A,I0)')       "INIT) -> SCAN : ", i, " NUMBER OF QM DATA POINTS : ", qmdata(i)%numdata
        WRITE(LOGUNIT,'(A,I0,A,I0)') "INIT) -> SCAN : ", i," NUMBER OF QM DATA POINTS : ", qmdata(i)%numdata
    END DO
    !*********************************************************************************************!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!********************************* CONFIGURATIONS RELATED **********************************!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! READ INFORMATION ABOUT CONFIGURATIONS FROM README FILE
    tempfile = ALLSPACE_FILE
    tempfile = TRIM(MDHEADDIR)//'/'//TRIM(conf_info_file)

    WRITE(LOGUNIT,'(A,A)') "INIT) READING CONFIGURATIONS INFO FILE : ", TRIM(tempfile)
    WRITE(*,'(A,A)')       "INIT) READING CONFIGURATIONS INFO FILE : ", TRIM(tempfile)
    CALL READCONFS(tempfile,myconfs)
    CALL SORTCONFS(myconfs)

    DO i = 1, NUMSCANS
        WRITE(*,'(A,I0,A,I0)')       "INIT) -> SCAN : ", i, " NUMBER OF CONFIGURATIONS : ", myconfs(i)%numconfs
        WRITE(LOGUNIT,'(A,I0,A,I0)') "INIT) -> SCAN : ", i," NUMBER OF CONFIGURATIONS : ", myconfs(i)%numconfs
    END DO
    !*********************************************************************************************!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! A CHECK TO ENSURE THAT WE ARE COMPARING THE QM AND MD DATA AT THE SAME DISTANCE POINTS
    CALL CHECK_QM_CONF_DIST(qmdata,myconfs)
    WRITE(*,'(A)')       "INIT) CONSISTANCY OF QM AND MD CONFS CHECKED.. ALL OK!"
    WRITE(LOGUNIT,'(A)') "INIT) CONSISTANCY OF QM AND MD CONFS CHECKED.. ALL OK!"
    !*********************************************************************************************!
 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!*********************************** MD DATA RELATED ***************************************!!    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! INITIALIZE MD DATA STRUCTURE AND WEIGHTS (SAME TYPE AS QM DATA)    
    WRITE(*,'(A)')       "INIT) CREATING MD DATA STRUCTURES.. "
    WRITE(LOGUNIT,'(A)') "INIT) CREATING MD DATA STRUCTURES.. "

    DO i = 1, NUMSCANS
        mdfirst(i)%numdata = qmdata(i)%numdata
        mdnow(i)%numdata   = qmdata(i)%numdata ;  mdnext(i)%numdata = qmdata(i)%numdata

        ALLOCATE(mdfirst(i)%qdist(mdfirst(i)%numdata),mdfirst(i)%qenergy(mdfirst(i)%numdata))
        ALLOCATE(mdnow(i)%qdist(mdnow(i)%numdata),mdnow(i)%qenergy(mdnow(i)%numdata))
        ALLOCATE(mdnext(i)%qdist(mdnext(i)%numdata),mdnext(i)%qenergy(mdnext(i)%numdata))
        ALLOCATE(mdfirst(i)%qweight(mdfirst(i)%numdata))
        ALLOCATE(mdnow(i)%qweight(mdnow(i)%numdata))
        ALLOCATE(mdnext(i)%qweight(mdnext(i)%numdata))

        mdfirst(i)%qdist(:) = qmdata(i)%qdist(:)
        mdnow(i)%qdist(:) = qmdata(i)%qdist(:) ;  mdnext(i)%qdist(:) = qmdata(i)%qdist(:)
        mdfirst(i)%qweight(:) = qmdata(i)%qweight(:)
        mdnow(i)%qweight(:) = qmdata(i)%qweight(:) ;  mdnext(i)%qweight(:) = qmdata(i)%qweight(:)
    END DO
    !*********************************************************************************************!
    
    !!***********************************************************************************************************************!!
    !!************************************************* INITIALIZATION PART *************************************************!!
    !!***********************************************************************************************************************!!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!************************************************** MONTE CARLO LOOP ***************************************************!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    WRITE(*,'(A)')       "!******************** STARTING MC LOOP ********************!"
    WRITE(LOGUNIT,'(A)') "!******************** STARTING MC LOOP ********************!"
    WRITE(*,'(A)')       "MCLOOP) INITIALIZATION.. "
    WRITE(LOGUNIT,'(A)') "MCLOOP) INITIALIZATION.. "
    
    109    FORMAT(I8,2X,F14.6)
    WRITE(myformat,'(A,I0,A)') "(I8,2X,",NUMLJPAR,"(F12.6,2X,F12.6,2X))"
    WRITE(myformat2,'(A,I0,A)') "(I8,2X,I8,2X",NUMLJPAR,"(F12.6,2X,F12.6,2X))"

    mciter = 0  ; mcaccept = 0 ; mcreject = 0 ; ACCEPTED = .FALSE. ; IACC = 0

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!************************************* 0TH MC STEP *****************************************!!    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    WRITE(*,'(A)')       "MCLOOP) DOING 0TH STEP.. "
    WRITE(LOGUNIT,'(A)') "MCLOOP) DOING 0TH STEP.. "
    
    CALL RUN_MD(mciter,firstitppar,mdfirst,firstcumerr,IMC_EQUIL)
    mdnow = mdfirst
    currcumerr = firstcumerr
    mincumerr = firstcumerr
    minpar(:) = firstpar(:)

    ! COPY guesspar.gtp to eq-0 

    WRITE(*,'(A,F14.6)')       "MCLOOP) >>>>>>>>>> FIRST RMSE: ", firstcumerr
    WRITE(LOGUNIT,'(A,F14.6)') "MCLOOP) >>>>>>>>>> FIRST RMSE: ", firstcumerr

    CALL INIT_MC_LOGS(IMC_EQUIL)
    CALL GETPARARR(firstpar,pararr)

    WRITE(MCALLENER,109) mciter, firstcumerr
    WRITE(MCACCENER,109) mciter, firstcumerr
    WRITE(MCREJENER,109) mciter, firstcumerr

    WRITE(MCALLPAR,TRIM(myformat)) mciter, pararr(:)
    WRITE(MCACCPAR,TRIM(myformat2)) mcaccept, mciter, pararr(:)
    WRITE(MCREJPAR,TRIM(myformat2)) mcreject, mciter, pararr(:)
    !*********************************************************************************************!

    WRITE(*,'(A)')       "!**********************************************************!"
    WRITE(LOGUNIT,'(A)') "!**********************************************************!"

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!************************************ EQUILIBRATION ****************************************!!    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    DO mciter = 1, MCEQUIL
        WRITE(*,'(A,I10,A)')       "!!!!!!!!!!!!!!!!!!!! MC EQUIL LOOP ", mciter," !!!!!!!!!!!!!!!!!!!!"
        WRITE(LOGUNIT,'(A,I10,A)') "!!!!!!!!!!!!!!!!!!!! MC EQUIL LOOP ", mciter," !!!!!!!!!!!!!!!!!!!!"

        ! MAKE A MOVE
        CALL MCMOVE(currpar,nextpar,nextitppar)

        ! RUN_MD
        CALL RUN_MD(mciter,nextitppar,mdnext,nextcumerr,IMC_EQUIL)
        delcumerr = nextcumerr - currcumerr
        IF(nextcumerr .LT. mincumerr) THEN
            mincumerr = nextcumerr
            minpar(:) = nextpar(:)
            miniter(1) = 0
            miniter(2) = mciter
        END IF
        
        WRITE(*,'(A,F14.6)')       "MCLOOP) BETA       : ", mcbeta%binit
        WRITE(LOGUNIT,'(A,F14.6)') "MCLOOP) BETA       : ", mcbeta%binit
        WRITE(*,'(A,F14.6)')       "MCLOOP) RMSE       : ", nextcumerr
        WRITE(LOGUNIT,'(A,F14.6)') "MCLOOP) RMSE       : ", nextcumerr
        WRITE(*,'(A,F14.6)')       "MCLOOP) DEL-RMSE   : ", delcumerr
        WRITE(LOGUNIT,'(A,F14.6)') "MCLOOP) DEL-RMSE   : ", delcumerr
        WRITE(*,'(A,F14.6)')       "MCLOOP) DELRMSE(%) : ", delcumerr*100.D0/currcumerr
        WRITE(LOGUNIT,'(A,F14.6)') "MCLOOP) DELRMSE(%) : ", delcumerr*100.D0/currcumerr

        ! ACCEPT/REJECT
        IF(delcumerr .LT. 0.D0) THEN
            ! ACCEPT
            ACCEPTED = .TRUE. ; IACC = 1
            mcaccept = mcaccept + 1
            currpar(:) = nextpar(:) ; currcumerr = nextcumerr

            WRITE(*,'(A)')       "MCLOOP)       VVVVV  ACCEPTED!  VVVVV"
            WRITE(LOGUNIT,'(A)') "MCLOOP)       VVVVV  ACCEPTED!  VVVVV"
        ELSE
            CALL RANDOM_NUMBER(myrandnum(1:1))
            IF (myrandnum(1) .LT. EXP(-1.D0*mcbeta%binit*delcumerr)) THEN
                ! ACCEPT
                ACCEPTED = .TRUE. ; IACC = 1
                mcaccept = mcaccept + 1
                currpar(:) = nextpar(:) ; currcumerr = nextcumerr
                
                WRITE(*,'(A)')       "MCLOOP)       VVVVV  ACCEPTED!  VVVVV"
                WRITE(LOGUNIT,'(A)') "MCLOOP)       VVVVV  ACCEPTED!  VVVVV"
            ELSE
                ! REJECT
                ACCEPTED = .FALSE. ; IACC = 0
                mcreject = mcreject + 1
                WRITE(*,'(A)')       "MCLOOP)       XXXXX  REJECTED!  XXXXX"
                WRITE(LOGUNIT,'(A)') "MCLOOP)       XXXXX  REJECTED!  XXXXX"
            ENDIF
        END IF

        ! LOGS
        CALL GETPARARR(currpar,pararr)
        WRITE(MCALLENER,109) mciter, currcumerr
        WRITE(MCALLPAR,TRIM(myformat)) mciter, pararr(:)

        
        IF(ACCEPTED) THEN
            WRITE(MCACCENER,109) mcaccept, nextcumerr
            WRITE(MCACCPAR,TRIM(myformat2)) mcaccept, mciter, pararr(:)
        ELSE
            CALL GETPARARR(nextpar,pararr)       
            WRITE(MCREJENER,109) mcreject, nextcumerr
            WRITE(MCREJPAR,TRIM(myformat2)) mcreject, mciter, pararr(:)
        ENDIF
        WRITE(*,'(A)')       "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        WRITE(LOGUNIT,'(A)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    END DO

    WRITE(*,'(A)')       "EQUIL) MONTE CARLO EQUILIBRATION LOOP ENEDED WITH THE FOLLOWING STATISTICS : "
    WRITE(LOGUNIT,'(A)') "EQUIL) MONTE CARLO EQUILIBRATION LOOP ENEDED WITH THE FOLLOWING STATISTICS : "

    WRITE(*,'(A,I8)')   "EQUIL) TOTAL MC LOOPS    : ", MCEQUIL
    WRITE(*,'(A,I8)')   "EQUIL) ACCEPTED MC LOOPS : ", mcaccept
    WRITE(*,'(A,I8)')   "EQUIL) REJECTED MC LOOPS : ", mcreject
    WRITE(*,'(A,F7.2)') "EQUIL) ACCEPTANCE PERCENTAGE : ", REAL(mcaccept*100)/REAL(MCEQUIL)
 
    WRITE(LOGUNIT,'(A,I8)')   "EQUIL) TOTAL MC LOOPS    : ", MCEQUIL
    WRITE(LOGUNIT,'(A,I8)')   "EQUIL) ACCEPTED MC LOOPS : ", mcaccept
    WRITE(LOGUNIT,'(A,I8)')   "EQUIL) REJECTED MC LOOPS : ", mcreject
    WRITE(LOGUNIT,'(A,F7.2)') "EQUIL) ACCEPTANCE PERCENTAGE : ", REAL(mcaccept*100)/REAL(MCEQUIL)
    CLOSE(MCALLENER) ; CLOSE(MCACCENER) ; CLOSE(MCREJENER)
    CLOSE(MCALLPAR)  ; CLOSE(MCACCPAR)  ; CLOSE(MCREJPAR)
    !*********************************************************************************************!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!************************************** QUENCHING ******************************************!!    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! CURRENT PAR (currpar) IS THE LAST ACCEPTED PARAMTER OF THE QUENCHING RUN
    mciter = 0 ; mcaccept = 0 ; mcreject = 0 ; ACCEPTED = .FALSE. ; IACC = 0
    CALL INITSENS()

    IF(MIN_PROD_START) THEN
        currpar(:) = minpar(:)
        currcumerr = mincumerr
        WRITE(*,'(A,I0)')       "MCPROD) PRODUCTION WILL BE STARTING FROM MINIMUM PAR SET OF EQUILIBRATION: eq-",&
                                miniter(2)
        WRITE(*,'(A,F14.6)')    "MCPROD) START RMSE     : ", firstcumerr
        WRITE(*,'(A,F14.6)')    "MCPROD) BEST EQUIL RMSE: ", mincumerr

        WRITE(LOGUNIT,'(A,I0)') "MCPROD) PRODUCTION WILL BE STARTING FROM MINIMUM PAR SET OF EQUILIBRATION: eq-",&
                                miniter(2)
        WRITE(LOGUNIT,'(A,F14.6)') "MCPROD) START RMSE     : ", firstcumerr
        WRITE(LOGUNIT,'(A,F14.6)') "MCPROD) BEST EQUIL RMSE: ", mincumerr
    END IF

    CALL INIT_MC_LOGS(IMC_PROD)  ; CALL INIT_ACC_LOGS()
    CALL GETPARARR(currpar,pararr)
    CALL GETSENSARR(currpar,sensarr)
    
    WRITE(MCALLENER,109) mciter, currcumerr
    WRITE(MCACCENER,109) mciter, currcumerr
    WRITE(MCREJENER,109) mciter, currcumerr

    WRITE(UNIT_SENSLOG,TRIM(myformat)) mciter, sensarr(:)
    WRITE(MCALLPAR,TRIM(myformat)) mciter, pararr(:)
    WRITE(MCACCPAR,TRIM(myformat2)) mcaccept, mciter, pararr(:)
    WRITE(MCREJPAR,TRIM(myformat2)) mcreject, mciter, pararr(:)

    ! --------------------------------------------------- MC LOOP STARTS -------------------------------------------------!
    DO mciter = 1, MCPROD
        WRITE(*,'(A,I10,A)')       "!******************* MC LOOP ", mciter," *******************!"
        WRITE(LOGUNIT,'(A,I10,A)') "!******************* MC LOOP ", mciter," *******************!"

        ! UPDATE BETA
        CALL GETBETA(mciter)

        ! UPDATE SENS
        IF(MOD(mciter,CONTROL_ACC_EVERY) == 0) THEN
            ! COLLECT ACCEPTANCE RATIO INFORMATION
            curr_acc_ratio = REAL(curr_acc_iter,KIND=kind_real)/REAL(CONTROL_ACC_EVERY,KIND=kind_real)
            err_acc_ratio = curr_acc_ratio - ref_acc_ratio
            interr_acc_ratio = interr_acc_ratio + err_acc_ratio
            diferr_acc_ratio = (err_acc_ratio-preverr_acc_ratio)

            ! PRINT ACC-RAT STATISTICS
            WRITE(UNIT_ACCLOG,'(I6,4(2X,F13.6))') mciter, curr_acc_ratio, err_acc_ratio,&
                                                  interr_acc_ratio, diferr_acc_ratio

            WRITE(*,'(A)')       "MCLOOP) ACCEPTANCE RATIO STATISTICS : "
            WRITE(LOGUNIT,'(A)') "MCLOOP) ACCEPTANCE RATIO STATISTICS : "

            WRITE(*,'(A,T2,A,T14,A,T29,A,T44,A,T59,A)') "#", " ITER ", " ACC-RAT ", &
                                                    " DEL-ACC ", " INT-ERR ", " DIF-ERR "

            WRITE(*,'(A,T2,A,T14,A,T29,A,T44,A,T59,A)') "#", " ---- ", " ------- ", &
                                                    " ------- ", " ------- ", " ------- "

            WRITE(LOGUNIT,'(A,T2,A,T14,A,T29,A,T44,A,T59,A)') "#", " ITER ", " ACC-RAT ", &
                                                    " DEL-ACC ", " INT-ERR ", " DIF-ERR "

            WRITE(LOGUNIT,'(A,T2,A,T14,A,T29,A,T44,A,T59,A)') "#", " ---- ", " ------- ", &
                                                    " ------- ", " ------- ", " ------- "
            
            WRITE(*,'(I6,4(2X,F13.6))') mciter, curr_acc_ratio, err_acc_ratio,&
                                                  interr_acc_ratio, diferr_acc_ratio

            WRITE(LOGUNIT,'(I6,4(2X,F13.6))') mciter, curr_acc_ratio, err_acc_ratio,&
                                                  interr_acc_ratio, diferr_acc_ratio
            
            WRITE(*,'(A,T2,A,T14,A,T29,A,T44,A,T59,A)') "#", " ---- ", " ------- ", &
                                                    " ------- ", " ------- ", " ------- "

            WRITE(LOGUNIT,'(A,T2,A,T14,A,T29,A,T44,A,T59,A)') "#", " ---- ", " ------- ", &
                                                    " ------- ", " ------- ", " ------- "

            ! UPDATE SENS
            IF (CONTROL_ACC) CALL UPDATESENS()

            ! WRITE TO GTP FILE
            currgtpfilename = ALLSPACE_FILE
            WRITE(currgtpfilename,'(A,A,I0,A)') TRIM(SENSDIR),'ps-',mciter,'.gtp'
            CALL WRITEGTP(currgtpfilename,curritppar)

            CALL GETSENSARR(currpar,sensarr)
            WRITE(UNIT_SENSLOG,TRIM(myformat)) mciter, sensarr(:)

            ! RESET ACC-RAT VARIABLES
            prev_acc_ratio = curr_acc_ratio ; curr_acc_iter = 0 ; preverr_acc_ratio = err_acc_ratio
        END IF

        ! MAKE A MOVE
        CALL MCMOVE(currpar,nextpar,nextitppar)

        ! RUN_MD
        CALL RUN_MD(mciter,nextitppar,mdnext,nextcumerr,IMC_PROD)
        delcumerr = nextcumerr - currcumerr

        ! KEEPING TRACK OF MINIMUM ERROR PARAMETER
        IF(nextcumerr .LT. mincumerr) THEN
            mincumerr = nextcumerr
            minpar(:) = nextpar(:)
            miniter(1) = 1
            miniter(2) = mciter
        END IF

        ! PRINT ESSENTIAL INFO        
        WRITE(*,'(A,F14.6)')       "MCLOOP) BETA       : ", mcbeta%bcurr
        WRITE(LOGUNIT,'(A,F14.6)') "MCLOOP) BETA       : ", mcbeta%bcurr
        WRITE(*,'(A,F14.6)')       "MCLOOP) RMSE       : ", nextcumerr
        WRITE(LOGUNIT,'(A,F14.6)') "MCLOOP) RMSE       : ", nextcumerr
        WRITE(*,'(A,F14.6)')       "MCLOOP) DEL-RMSE   : ", delcumerr
        WRITE(LOGUNIT,'(A,F14.6)') "MCLOOP) DEL-RMSE   : ", delcumerr
        WRITE(*,'(A,F14.6)')       "MCLOOP) DELRMSE(%) : ", delcumerr*100.D0/currcumerr
        WRITE(LOGUNIT,'(A,F14.6)') "MCLOOP) DELRMSE(%) : ", delcumerr*100.D0/currcumerr

        ! ACCEPT/REJECT
        IF(delcumerr .LT. 0.D0) THEN
            ! ACCEPT
            ACCEPTED = .TRUE. ; IACC = 1
            mcaccept = mcaccept + 1 ; curr_acc_iter = curr_acc_iter + 1
            currpar(:) = nextpar(:) ; currcumerr = nextcumerr

            WRITE(*,'(A)')       "MCLOOP)       VVVVV  ACCEPTED!  VVVVV"
            WRITE(LOGUNIT,'(A)') "MCLOOP)       VVVVV  ACCEPTED!  VVVVV"
        ELSE
            CALL RANDOM_NUMBER(myrandnum(1:1))
            IF (myrandnum(1) .LT. EXP(-1.D0*mcbeta%bcurr*delcumerr)) THEN
                ! ACCEPT
                ACCEPTED = .TRUE. ; IACC = 1
                mcaccept = mcaccept + 1 ; curr_acc_iter = curr_acc_iter + 1 
                currpar(:) = nextpar(:) ; currcumerr = nextcumerr
                
                WRITE(*,'(A)')       "MCLOOP)       VVVVV  ACCEPTED!  VVVVV"
                WRITE(LOGUNIT,'(A)') "MCLOOP)       VVVVV  ACCEPTED!  VVVVV"
            ELSE
                ! REJECT
                ACCEPTED = .FALSE. ; IACC = 0
                mcreject = mcreject + 1
                WRITE(*,'(A)')       "MCLOOP)       XXXXX  REJECTED!  XXXXX"
                WRITE(LOGUNIT,'(A)') "MCLOOP)       XXXXX  REJECTED!  XXXXX"
            ENDIF
        END IF

        ! LOGS
        CALL GETPARARR(currpar,pararr)
        WRITE(MCALLENER,109) mciter, currcumerr
        WRITE(MCALLPAR,TRIM(myformat)) mciter, pararr(:)
        IF(ACCEPTED) THEN
            WRITE(MCACCENER,109) mcaccept, nextcumerr
            WRITE(MCACCPAR,TRIM(myformat2)) mcaccept, mciter, pararr(:)
        ELSE
            CALL GETPARARR(nextpar,pararr)       
            WRITE(MCREJENER,109) mcreject, nextcumerr
            WRITE(MCREJPAR,TRIM(myformat2)) mcreject, mciter, pararr(:)
        ENDIF

        WRITE(*,'(A)')       "!**********************************************************!"
        WRITE(LOGUNIT,'(A)') "!**********************************************************!"
    END DO
    ! --------------------------------------------------- MC LOOP ENDS -------------------------------------------------!
    
    WRITE(*,'(A)')       "MCPROD) MONTE CARLO LOOP ENEDED WITH THE FOLLOWING STATISTICS : "
    WRITE(LOGUNIT,'(A)') "MCPROD) MONTE CARLO LOOP ENEDED WITH THE FOLLOWING STATISTICS : "

    WRITE(*,'(A,I8)')   "MCPROD) TOTAL MC LOOPS    : ", MCPROD
    WRITE(*,'(A,I8)')   "MCPROD) ACCEPTED MC LOOPS : ", mcaccept
    WRITE(*,'(A,I8)')   "MCPROD) REJECTED MC LOOPS : ", mcreject
    WRITE(*,'(A,F7.2)') "MCPROD) ACCEPTANCE PERCENTAGE : ", REAL(mcaccept*100)/REAL(MCPROD)
 
    WRITE(LOGUNIT,'(A,I8)')   "MCPROD) TOTAL MC LOOPS    : ", MCPROD
    WRITE(LOGUNIT,'(A,I8)')   "MCPROD) ACCEPTED MC LOOPS : ", mcaccept
    WRITE(LOGUNIT,'(A,I8)')   "MCPROD) REJECTED MC LOOPS : ", mcreject
    WRITE(LOGUNIT,'(A,F7.2)') "MCPROD) ACCEPTANCE PERCENTAGE : ", REAL(mcaccept*100)/REAL(MCPROD)
    CLOSE(MCALLENER) ; CLOSE(MCACCENER) ; CLOSE(MCREJENER)
    CLOSE(MCALLPAR)  ; CLOSE(MCACCPAR)  ; CLOSE(MCREJPAR)
    !*********************************************************************************************!

    WRITE(*,*)
    WRITE(LOGUNIT,*)

    WRITE(*,'(A)')       "!!!!!!........................................................................................!!!!!!"
    WRITE(LOGUNIT,'(A)') "!!!!!!........................................................................................!!!!!!"

    WRITE(*,'(A)')       "!!!!!!........................................................................................!!!!!!"
    WRITE(LOGUNIT,'(A)') "!!!!!!........................................................................................!!!!!!"

    WRITE(*,'(A,T26,A,T99,A)')       "!!","TUNELJ) TUNED LJ PARAMETER DETAILS : ","!!"
    WRITE(LOGUNIT,'(A,T26,A,T99,A)') "!!","TUNELJ) TUNED LJ PARAMETER DETAILS : ","!!"

    WRITE(*,'(A,T26,A,F14.6,T99,A)')       "!!","TUNELJ) FIRST RMSE         : ", firstcumerr,"!!"
    WRITE(*,'(A,T26,A,F14.6,T99,A)')       "!!","TUNELJ) MIN   RMSE         : ", mincumerr,"!!"
    WRITE(*,'(A,T26,A,F14.6,T99,A)')       "!!","TUNELJ) DEL MIN RMSE       : ", mincumerr-firstcumerr,"!!"
    WRITE(*,'(A,T26,A,F14.6,T99,A)')       "!!","TUNELJ) DEL MIN  RMSE(%)   : ", (mincumerr-firstcumerr)*100.D0/firstcumerr,"!!"
    WRITE(*,'(A,T26,A,F14.6,T99,A)')       "!!","TUNELJ) LAST  RMSE         : ", currcumerr,"!!"
    WRITE(*,'(A,T26,A,F14.6,T99,A)')       "!!","TUNELJ) DEL LAST RMSE      : ", currcumerr-firstcumerr,"!!"
    WRITE(*,'(A,T26,A,F14.6,T99,A)')       "!!","TUNELJ) DEL LAST  RMSE(%)  : ", (currcumerr-firstcumerr)*100.D0/firstcumerr,"!!"

    WRITE(LOGUNIT,'(A,T26,A,F14.6,T99,A)') "!!","TUNELJ) FIRST RMSE         : ", firstcumerr,"!!"
    WRITE(LOGUNIT,'(A,T26,A,F14.6,T99,A)') "!!","TUNELJ) MIN   RMSE         : ", mincumerr,"!!"
    WRITE(LOGUNIT,'(A,T26,A,F14.6,T99,A)') "!!","TUNELJ) DEL MIN RMSE       : ", mincumerr-firstcumerr,"!!"
    WRITE(LOGUNIT,'(A,T26,A,F14.6,T99,A)') "!!","TUNELJ) DEL MIN  RMSE(%)   : ", (mincumerr-firstcumerr)*100.D0/firstcumerr,"!!"
    WRITE(LOGUNIT,'(A,T26,A,F14.6,T99,A)') "!!","TUNELJ) LAST  RMSE         : ", currcumerr,"!!"
    WRITE(LOGUNIT,'(A,T26,A,F14.6,T99,A)') "!!","TUNELJ) DEL LAST RMSE      : ", currcumerr-firstcumerr,"!!"
    WRITE(LOGUNIT,'(A,T26,A,F14.6,T99,A)') "!!","TUNELJ) DEL LAST  RMSE(%)  : ", (currcumerr-firstcumerr)*100.D0/firstcumerr,"!!"

    IF(miniter(1) == 0) THEN
        WRITE(*,'(A,T26,A,I0,T99,A)')      "!!","TUNELJ) MINIMUM RMSE FRAME : eq-", miniter(2),"!!"
        WRITE(LOGUNIT,'(A,T26,A,I0,T99,A)')    "!!","TUNELJ) MINIMUM RMSE FRAME : eq-", miniter(2),"!!"
    ELSE
        WRITE(*,'(A,T26,A,I0,T99,A)')    "!!","TUNELJ) MINIMUM RMSE FRAME : ps-", miniter(2),"!!"
        WRITE(LOGUNIT,'(A,T26,A,I0,T99,A)')    "!!","TUNELJ) MINIMUM RMSE FRAME : ps-", miniter(2),"!!"
    END IF
    WRITE(*,'(A,T99,A)') "!!", "!!"  
    WRITE(LOGUNIT,'(A,T99,A)') "!!", "!!"    

    
    WRITE(LOGUNIT,*)

    WRITE(*,'(A,A,T15,A,T30,A,T43,A,T57,A,T71,A,T85,A,T99,A))') "!!"," --> PAR "," INIT-SIG "," MIN-SIG "," DEL-SIG(%) ", &
                                            " INIT-EPS "," MIN-EPS "," DEL-EPS(%) ","!!"
    WRITE(*,'(A,A,T15,A,T30,A,T43,A,T57,A,T71,A,T85,A,T99,A))') "!!"," --> --- "," -------- "," ------- "," ---------- ", &
                                            " -------- "," ------- "," ---------- ","!!"
    WRITE(LOGUNIT,'(A,A,T15,A,T30,A,T43,A,T57,A,T71,A,T85,A,T99,A))') "!!","--> PAR "," INIT-SIG "," MIN-SIG "," DEL-SIG(%) ", &
                                            " INIT-EPS "," MIN-EPS "," DEL-EPS(%) ","!!"
    WRITE(LOGUNIT,'(A,A,T15,A,T30,A,T43,A,T57,A,T71,A,T85,A,T99,A))') "!!"," --> --- "," -------- "," ------- "," ---------- ", &
                                            " -------- "," ------- "," ---------- ","!!"
    DO i = 1, NUMLJPAR
        WRITE(*,'(A,A,A,2X,6(F12.6,2X),T99,A)') "!!"," --> ", TRIM(firstpar(i)%parname), &
                                firstpar(i)%par(1), minpar(i)%par(1), &
                                (minpar(i)%par(1) - firstpar(i)%par(1))*100.D0/firstpar(i)%par(1), &
                                firstpar(i)%par(2), minpar(i)%par(2), &
                                (minpar(i)%par(2) - firstpar(i)%par(2))*100.D0/firstpar(i)%par(2),"!!"
        WRITE(LOGUNIT,'(A,A,A,2X,6(F12.6,2X),T99,A)') "!!"," --> ", TRIM(firstpar(i)%parname), &
                                firstpar(i)%par(1), minpar(i)%par(1), &
                                (minpar(i)%par(1) - firstpar(i)%par(1))*100.D0/firstpar(i)%par(1), &
                                firstpar(i)%par(2), minpar(i)%par(2), &
                                (minpar(i)%par(2) - firstpar(i)%par(2))*100.D0/firstpar(i)%par(2),"!!"
    END DO
    WRITE(*,'(A,T99,A)') "!!", "!!"    
    WRITE(LOGUNIT,'(A,T99,A)') "!!", "!!"    

    WRITE(*,'(A,A,T15,A,T30,A,T43,A,T57,A,T71,A,T85,A,T99,A))') "!!"," --> PAR "," INIT-SIG "," LAST-SIG "," DEL-SIG(%) ", &
                                            " INIT-EPS "," LAST-EPS "," DEL-EPS(%) ","!!"
    WRITE(LOGUNIT,'(A,A,T15,A,T30,A,T43,A,T57,A,T71,A,T85,A,T99,A))') "!!"," --> PAR "," INIT-SIG "," LAST-SIG "," DEL-SIG(%) ", &
                                            " INIT-EPS "," LAST-EPS "," DEL-EPS(%) ","!!"
    WRITE(*,'(A,A,T15,A,T30,A,T43,A,T57,A,T71,A,T85,A,T99,A))') "!!"," --> --- "," -------- "," -------- "," ---------- ", &
                                            " -------- "," -------- "," ---------- ","!!"
    WRITE(LOGUNIT,'(A,A,T15,A,T30,A,T43,A,T57,A,T71,A,T85,A,T99,A))') "!!"," --> --- "," -------- "," -------- "," ---------- ", &
                                            " -------- "," -------- "," ---------- ","!!"
    DO i = 1, NUMLJPAR
        WRITE(*,'(A,A,A,2X,6(F12.6,2X),T99,A)') "!!"," --> ", TRIM(firstpar(i)%parname), &
                                firstpar(i)%par(1), currpar(i)%par(1), &
                                (currpar(i)%par(1) - firstpar(i)%par(1))*100.D0/firstpar(i)%par(1), &
                                firstpar(i)%par(2), currpar(i)%par(2), &
                                (currpar(i)%par(2) - firstpar(i)%par(2))*100.D0/firstpar(i)%par(2),"!!"
        WRITE(LOGUNIT,'(A,A,A,2X,6(F12.6,2X),T99,A)') "!!"," --> ", TRIM(firstpar(i)%parname), &
                                firstpar(i)%par(1), currpar(i)%par(1), &
                                (currpar(i)%par(1) - firstpar(i)%par(1))*100.D0/firstpar(i)%par(1), &
                                firstpar(i)%par(2), currpar(i)%par(2), &
                                (currpar(i)%par(2) - firstpar(i)%par(2))*100.D0/firstpar(i)%par(2),"!!"
    END DO

    WRITE(*,'(A)')       "!!!!!!........................................................................................!!!!!!"
    WRITE(LOGUNIT,'(A)') "!!!!!!........................................................................................!!!!!!"
    WRITE(*,'(A)')       "!!!!!!........................................................................................!!!!!!"
    WRITE(LOGUNIT,'(A)') "!!!!!!........................................................................................!!!!!!"

    !!***********************************************************************************************************************!!
    !!************************************************** MONTE CARLO LOOP ***************************************************!!
    !!***********************************************************************************************************************!!

    WRITE(*,'(T40,A)')       " >>> ALL THE BEST :) <<< "
    WRITE(LOGUNIT,'(T40,A)') " >>> ALL THE BEST :) <<< "

    CONTAINS

    SUBROUTINE PARSEINP(filename)
        CHARACTER(LEN=len_filename), INTENT(IN) :: filename
        CHARACTER(LEN=len_filename) :: ctemp
        INTEGER :: ios, INPUNIT
        CHARACTER(LEN=len_record) :: emsg


        OPEN(NEWUNIT=INPUNIT,FILE=TRIM(filename),IOSTAT=ios,STATUS='OLD',ACTION='READ')
        emsg = "PARSEINP) ERROR: CANNOT OPEN INPUT FILE : "//TRIM(filename)
        CALL CHECKIOS(ios,emsg)

        emsg = "PARSEINP) ERROR: WHILE READING INPUT FILE : "//TRIM(filename)

        READ(INPUNIT,*,IOSTAT=ios)
        READ(INPUNIT,*,IOSTAT=ios) NUMLJPAR
        CALL CHECKIOS(ios,emsg)

        READ(INPUNIT,*,IOSTAT=ios)
        READ(INPUNIT,*,IOSTAT=ios) NUMSCANS
        CALL CHECKIOS(ios,emsg)

        READ(INPUNIT,*,IOSTAT=ios)
        READ(INPUNIT,'(A)',IOSTAT=ios) HEADDIR
        CALL CHECKIOS(ios,emsg)
        ctemp = ADJUSTL(HEADDIR)
        HEADDIR = ALLSPACE_FILE
        HEADDIR = ctemp

        READ(INPUNIT,*,IOSTAT=ios)
        READ(INPUNIT,'(A)',IOSTAT=ios) MDHEADDIR
        CALL CHECKIOS(ios,emsg)
        ctemp = ADJUSTL(MDHEADDIR)
        MDHEADDIR = ALLSPACE_FILE
        MDHEADDIR = ctemp
        
        READ(INPUNIT,*,IOSTAT=ios)
        READ(INPUNIT,'(A)',IOSTAT=ios) QMHEADDIR
        CALL CHECKIOS(ios,emsg)
        ctemp = ADJUSTL(QMHEADDIR)
        QMHEADDIR = ALLSPACE_FILE
        QMHEADDIR = ctemp
        
        READ(INPUNIT,*,IOSTAT=ios)
        READ(INPUNIT,'(A)',IOSTAT=ios) qm_info_file
        CALL CHECKIOS(ios,emsg)
        ctemp = ADJUSTL(qm_info_file)
        qm_info_file = ALLSPACE_FILE
        qm_info_file = ctemp

        READ(INPUNIT,*,IOSTAT=ios)
        READ(INPUNIT,'(A)',IOSTAT=ios) conf_info_file
        CALL CHECKIOS(ios,emsg)
        ctemp = ADJUSTL(conf_info_file)
        conf_info_file = ALLSPACE_FILE
        conf_info_file = ctemp
        
        READ(INPUNIT,*,IOSTAT=ios)
        READ(INPUNIT,'(A)',IOSTAT=ios) mdpfilename
        CALL CHECKIOS(ios,emsg)
        ctemp = ADJUSTL(mdpfilename)
        mdpfilename = ALLSPACE_FILE
        mdpfilename = ctemp
        
        READ(INPUNIT,*,IOSTAT=ios)
        READ(INPUNIT,'(A)',IOSTAT=ios) gtpfilename
        CALL CHECKIOS(ios,emsg)
        ctemp = ADJUSTL(gtpfilename)
        gtpfilename = ALLSPACE_FILE
        gtpfilename = ctemp
        
        
        READ(INPUNIT,*,IOSTAT=ios)
        READ(INPUNIT,'(A)',IOSTAT=ios) GROMACSEXE
        CALL CHECKIOS(ios,emsg)
        ctemp = ADJUSTL(GROMACSEXE)
        GROMACSEXE = ALLSPACE_FILE
        GROMACSEXE = ctemp

        READ(INPUNIT,*,IOSTAT=ios)
        READ(INPUNIT,*,IOSTAT=ios) MCEQUIL
        CALL CHECKIOS(ios,emsg)

        READ(INPUNIT,*,IOSTAT=ios)
        READ(INPUNIT,*,IOSTAT=ios) MCPROD
        CALL CHECKIOS(ios,emsg)

        READ(INPUNIT,*,IOSTAT=ios)
        READ(INPUNIT,*,IOSTAT=ios) mcbeta%numpoints
        CALL CHECKIOS(ios,emsg)

        READ(INPUNIT,*,IOSTAT=ios)
        READ(INPUNIT,*,IOSTAT=ios) mcbeta%binit
        CALL CHECKIOS(ios,emsg)

        ALLOCATE(mcbeta%biter(0:mcbeta%numpoints+1),mcbeta%beta(0:mcbeta%numpoints+1))

        READ(INPUNIT,*,IOSTAT=ios)
        READ(INPUNIT,*,IOSTAT=ios) mcbeta%beta(1:mcbeta%numpoints)
        CALL CHECKIOS(ios,emsg)

        READ(INPUNIT,*,IOSTAT=ios)
        READ(INPUNIT,*,IOSTAT=ios) mcbeta%biter(1:mcbeta%numpoints)
        CALL CHECKIOS(ios,emsg)

        READ(INPUNIT,*,IOSTAT=ios)
        READ(INPUNIT,*,IOSTAT=ios) MIN_PROD_START
        emsg = ALLSPACE_REC 
        emsg = "PARSEINP) ERROR: COULD NOT READ MINIMUM PRODUCTION RESTART FLAG. GIVE .TRUE. OR .FALSE."
        CALL CHECKIOS(ios,emsg)

        READ(INPUNIT,*,IOSTAT=ios)
        READ(INPUNIT,*,IOSTAT=ios) CONTROL_ACC
        emsg = ALLSPACE_REC 
        emsg = "PARSEINP) ERROR: COULD NOT READ CONTROL ACCEPTANCE RATIO FLAG. GIVE .TRUE. OR .FALSE."
        CALL CHECKIOS(ios,emsg)

        READ(INPUNIT,*,IOSTAT=ios)
        READ(INPUNIT,*,IOSTAT=ios) CONTROL_ACC_EVERY
        emsg = ALLSPACE_REC 
        emsg = "PARSEINP) ERROR: WHILE READING INPUT FILE : "//TRIM(filename)
        CALL CHECKIOS(ios,emsg)

        READ(INPUNIT,*,IOSTAT=ios)
        READ(INPUNIT,*,IOSTAT=ios) ref_acc_ratio
        CALL CHECKIOS(ios,emsg)
        ! CONVERT PERCENT TO RATIO
        ref_acc_ratio = ref_acc_ratio*0.01D0

        READ(INPUNIT,*,IOSTAT=ios)
        READ(INPUNIT,*,IOSTAT=ios) kp, ki, kd, sens_scale
        CALL CHECKIOS(ios,emsg)

        READ(INPUNIT,*,IOSTAT=ios)
        READ(INPUNIT,*,IOSTAT=ios) MINSENS_SIG, MINSENS_EPS
        CALL CHECKIOS(ios,emsg)

        CLOSE(INPUNIT)
    END SUBROUTINE PARSEINP

    SUBROUTINE PRINT_LOG_PARAMS()
        WRITE(*,'(A)')    "INIT) HEADDIR                    : "//TRIM(HEADDIR)
        WRITE(*,'(A)')    "INIT) RUNDIR                     : "//TRIM(RUNDIR)
        WRITE(*,'(A)')    "INIT) MD HEADDIR                 : "//TRIM(MDHEADDIR)
        WRITE(*,'(A)')    "INIT) QM HEADDIR                 : "//TRIM(QMHEADDIR)
        WRITE(*,'(A)')    "INIT) QM SCANS FILENAME          : "//TRIM(qm_info_file)
        WRITE(*,'(A)')    "INIT) CONF DETAILS FILENAME      : "//TRIM(conf_info_file)
        WRITE(*,'(A)')    "INIT) MDP FILENAME               : "//TRIM(mdpfilename)
        WRITE(*,'(A)')    "INIT) GTP FILENAME               : "//TRIM(gtpfilename)
        WRITE(*,'(A,I0)') "INIT) NUMBER OF LJ PARAMTERS     : ", NUMLJPAR
        WRITE(*,'(A,I0)') "INIT) NUMBER OF SCANS            : ", NUMSCANS
        WRITE(*,'(A,I0)') "INIT) NUMBER OF ANNEALING POINTS : ", mcbeta%numpoints
        WRITE(*,'(A,I0)') "INIT) NUMBER OF MC EQUILIBRATION : ", MCEQUIL
        WRITE(*,'(A,I0)') "INIT) NUMBER OF MC PRODUCTION    : ", MCPROD
        WRITE(*,'(A,L1)') "INIT) STARTING PRODUCTION FOR MIN: ", MIN_PROD_START
    
        WRITE(LOGUNIT,'(A)')    "INIT) HEADDIR                    : "//TRIM(HEADDIR)
        WRITE(LOGUNIT,'(A)')    "INIT) RUNDIR                     : "//TRIM(RUNDIR)
        WRITE(LOGUNIT,'(A)')    "INIT) MD HEADDIR                 : "//TRIM(MDHEADDIR)
        WRITE(LOGUNIT,'(A)')    "INIT) QM HEADDIR                 : "//TRIM(QMHEADDIR)
        WRITE(LOGUNIT,'(A)')    "INIT) QM SCANS FILENAME          : "//TRIM(qm_info_file)
        WRITE(LOGUNIT,'(A)')    "INIT) CONF DETAILS FILENAME      : "//TRIM(conf_info_file)
        WRITE(LOGUNIT,'(A)')    "INIT) MDP FILENAME               : "//TRIM(mdpfilename)
        WRITE(LOGUNIT,'(A)')    "INIT) GTP FILENAME               : "//TRIM(gtpfilename)
        WRITE(LOGUNIT,'(A,I0)') "INIT) NUMBER OF LJ PARAMTERS     : ", NUMLJPAR
        WRITE(LOGUNIT,'(A,I0)') "INIT) NUMBER OF SCANS            : ", NUMSCANS
        WRITE(LOGUNIT,'(A,I0)') "INIT) NUMBER OF ANNEALING POINTS : ", mcbeta%numpoints
        WRITE(LOGUNIT,'(A,I0)') "INIT) NUMBER OF MC EQUILIBRATION : ", MCEQUIL
        WRITE(LOGUNIT,'(A,I0)') "INIT) NUMBER OF MC PRODUCTION    : ", MCPROD
        WRITE(LOGUNIT,'(A,L1)') "INIT) STARTING PRODUCTION FOR MIN: ", MIN_PROD_START
    END SUBROUTINE PRINT_LOG_PARAMS

    SUBROUTINE READQMDATA(filename,myqmdata)
        CHARACTER(LEN=len_filename), INTENT(INOUT) :: filename
        TYPE(qmdatatype), INTENT(INOUT)         :: myqmdata(:)
        INTEGER :: i, ios, itemp
        INTEGER :: QMUNIT, QMSCAN
        CHARACTER(LEN=len_filename) :: tempfilename 
        CHARACTER(LEN=len_record)   :: commentline
        CHARACTER(LEN=len_scanname) :: tempscanname

        OPEN(NEWUNIT=QMUNIT,FILE=TRIM(filename),IOSTAT=ios,STATUS='OLD',ACTION='READ')
        IF(ios .NE. 0) THEN
            WRITE(*,'(A,A)') "ERROR) CANNOT OPEN QM INFO FILE : ", TRIM(filename)
            STOP
        END IF

        READ(QMUNIT,'(A)',IOSTAT=ios) commentline
        WRITE(*,'(A,A)') "READQMDATA) COMMENT ON INFO FILE : ",TRIM(commentline)
        WRITE(LOGUNIT,'(A,A)') "READQMDATA) COMMENT ON INFO FILE : ",TRIM(commentline)
        
        ! READ THE NUMBER OF SCANS FROM THE SCANS INFO FILE AND COMPARE IT WITH JOB CONTROL INPUT
        READ(QMUNIT,*,IOSTAT=ios) itemp
        IF(itemp .NE. NUMSCANS) THEN
            WRITE(*,'(A)') " READQMDATA) ERROR : NUMSCANS IN QM INFO FILE DOESN'T MATCH JOB INPUT FILE'S VALUE " 
            WRITE(*,'(A,I0,A,I0)') " READQMDATA) QM INFO NUMSCANS : ", itemp,"JOB CONTROL NUMSCANS : ", NUMSCANS 
            WRITE(LOGUNIT,'(A)') " READQMDATA) ERROR : NUMSCANS IN QM INFO FILE DOESN'T MATCH JOB INPUT FILE'S VALUE " 
            WRITE(LOGUNIT,'(A,I0,A,I0)') " READQMDATA) QM INFO NUMSCANS : ", itemp,"JOB CONTROL NUMSCANS : ", NUMSCANS 
            STOP
        END IF

        ! READ SCANS INFO FILE FOR INFO(SCAN NAME || QM ENERGY FILE NAME) ABOUT EACH SCAN
        DO i = 1, NUMSCANS
            READ(QMUNIT,*,IOSTAT=ios) qmdata(i)%qscanname, qmdata(i)%qfilename
            IF(ios .NE. 0) THEN
                WRITE(*,'(A,A)') "ERROR) WHILE READING QM INFO FILE : ", TRIM(filename)
                STOP
            END IF

            ! ADJUST THE SCAN NAME AND FILE NAME
            tempscanname = ALLSPACE_SCAN
            tempscanname = TRIM(ADJUSTL(qmdata(i)%qscanname))
            qmdata(i)%qscanname = ALLSPACE_SCAN
            qmdata(i)%qscanname = TRIM(tempscanname)

            tempfilename = ALLSPACE_FILE
            tempfilename = TRIM(ADJUSTL(qmdata(i)%qfilename))
            qmdata(i)%qfilename = ALLSPACE_FILE
            qmdata(i)%qfilename = TRIM(tempfilename)

            tempfilename = ALLSPACE_FILE
            tempfilename = TRIM(QMHEADDIR)//'/'//TRIM(qmdata(i)%qfilename)

            ! COPY THE QM DATA INTO THE RUN DIRECTORY -- SERVES AS A LOG
            CALL FCOPY(tempfilename,DATADIR)

            ! OPEN THE QM ENERGY FILE CORRESPONDING TO A CERTAIN SCAN
            OPEN(NEWUNIT=QMSCAN,FILE=TRIM(tempfilename),IOSTAT=ios,STATUS='OLD',ACTION='READ')
            IF(ios .NE. 0) THEN
                WRITE(*,'(A,A)') "ERROR) CANNOT OPEN QM DATA FILE : ", TRIM(tempfilename)
                STOP
            END IF

            WRITE(*,'(A,A)') " READQMDATA) READING FROM QM DATA FILE : ",TRIM(tempfilename)
            WRITE(LOGUNIT,'(A,A)') " READQMDATA) READING FROM QM DATA FILE : ",TRIM(tempfilename)

            ! READ AND ALLOCATE THE NUMBER OF CONFIGURATIONS IN THIS SCAN            
            READ(QMSCAN,*,IOSTAT=ios) 
            READ(QMSCAN,*,IOSTAT=ios) qmdata(i)%numdata
            ALLOCATE(qmdata(i)%qdist(qmdata(i)%numdata),qmdata(i)%qenergy(qmdata(i)%numdata))
            ALLOCATE(qmdata(i)%qweight(qmdata(i)%numdata))

            ! READ DISTANCE, ENERGY AND THE CORRESPONDING WEIGHTS
            DO j = 1, qmdata(i)%numdata
                READ(QMSCAN,*,IOSTAT=ios) qmdata(i)%qdist(j), qmdata(i)%qenergy(j), qmdata(i)%qweight(j)
                IF(ios .NE. 0) THEN
                    WRITE(*,'(A,A)') "ERROR) WHILE READING QM FILE : ", TRIM(tempfilename)
                    STOP
                END IF
            END DO
            CLOSE(QMSCAN)    
        END DO 
        CLOSE(QMUNIT)
    END SUBROUTINE READQMDATA

    SUBROUTINE SORTQMDATA(myqmdata)
        TYPE(qmdatatype), INTENT(INOUT) :: myqmdata(:)
        INTEGER :: isc, ios, i, j, imintemp
        REAL(KIND=kind_real) :: rmintemp
        REAL(KIND=kind_real) :: qm_einf
        CHARACTER(LEN=len_record) :: emsg

        DO isc = 1, NUMSCANS
            IF( .NOT. ALLOCATED(qmdata(isc)%qdist)) STOP "INTERNAL ERROR"
            IF( .NOT. ALLOCATED(qmdata(isc)%qenergy)) STOP "INTERNAL ERROR"
            DO i = 1, qmdata(isc)%numdata-1
                rmintemp = qmdata(isc)%qdist(i) ; imintemp = i
                DO j = i+1, qmdata(isc)%numdata
                    IF(qmdata(isc)%qdist(j) .LT. rmintemp) THEN
                        rmintemp = qmdata(isc)%qdist(j)
                        imintemp = j
                    END IF
                END DO
                IF(imintemp .NE. i) THEN
                    rmintemp = qmdata(isc)%qdist(imintemp)
                    qmdata(isc)%qdist(imintemp) = qmdata(isc)%qdist(i)
                    qmdata(isc)%qdist(i) = rmintemp
                    rmintemp = qmdata(isc)%qenergy(imintemp)
                    qmdata(isc)%qenergy(imintemp) = qmdata(isc)%qenergy(i)
                    qmdata(isc)%qenergy(i) = rmintemp
                    rmintemp = qmdata(isc)%qweight(imintemp)
                    qmdata(isc)%qweight(imintemp) = qmdata(isc)%qweight(i)
                    qmdata(isc)%qweight(i) = rmintemp
                END IF
            END DO

            ! SHIFT QM ENERGY TO ZERO IF NOT DONE ALREADY
            qm_einf = qmdata(isc)%qenergy(qmdata(isc)%numdata)
            qmdata(isc)%qenergy(:) = qmdata(isc)%qenergy(:) - qm_einf
        END DO
    END SUBROUTINE SORTQMDATA

    SUBROUTINE READCONFS(filename,myconfs)
        CHARACTER(LEN=len_filename), INTENT(INOUT) :: filename
        TYPE(conftype), INTENT(INOUT)           :: myconfs(:)
        INTEGER :: isc, i, itemp, ios, CONFUNIT, CONFSCAN
        CHARACTER(LEN=len_record)   :: commentline
        CHARACTER(LEN=len_scanname) :: tempscanname
        CHARACTER(LEN=len_dirname) :: tempdirname

        OPEN(NEWUNIT=CONFUNIT,FILE=TRIM(filename),IOSTAT=ios,STATUS='OLD',ACTION='READ')
        IF(ios .NE. 0) THEN
            WRITE(*,'(A,A)') "ERROR) CANNOT OPEN CONF INFO FILE : ", TRIM(filename)
            STOP
        END IF

        READ(CONFUNIT,'(A)',IOSTAT=ios) commentline
        WRITE(*,'(A,A)') "READCONFS) COMMENT ON INFO FILE : ",TRIM(commentline)
        WRITE(LOGUNIT,'(A,A)') "READCONFS) COMMENT ON INFO FILE : ",TRIM(commentline)

        READ(CONFUNIT,*,IOSTAT=ios) itemp
        IF(itemp .NE. NUMSCANS) THEN
            WRITE(*,'(A)') " READCONFS) ERROR : NUMSCANS IN MD INFO FILE DOESN'T MATCH JOB INPUT FILE'S VALUE " 
            WRITE(*,'(A,I0,A,I0)') " READCONFS) MD INFO NUMSCANS : ", itemp,"JOB CONTROL NUMSCANS : ", NUMSCANS 
            WRITE(LOGUNIT,'(A)') " READCONFS) ERROR : NUMSCANS IN MD INFO FILE DOESN'T MATCH JOB INPUT FILE'S VALUE " 
            WRITE(LOGUNIT,'(A,I0,A,I0)') " READCONFS) MD INFO NUMSCANS : ", itemp,"JOB CONTROL NUMSCANS : ", NUMSCANS 
            STOP
        END IF

        DO isc = 1, NUMSCANS
            ! READ THE CONFIGURATIONS RELATED INFO (SCAN NAME || SCAN DIRECTORY || README FILE)
            READ(CONFUNIT,*,IOSTAT=ios) myconfs(isc)%cscanname, myconfs(isc)%cdirname, myconfs(isc)%cfilename

            tempscanname = ALLSPACE_SCAN
            tempscanname = TRIM(ADJUSTL(myconfs(isc)%cscanname))
            myconfs(isc)%cscanname = ALLSPACE_SCAN
            myconfs(isc)%cscanname = TRIM(tempscanname)

            tempdirname = ALLSPACE_DIR
            tempdirname = TRIM(ADJUSTL(myconfs(isc)%cdirname))
            myconfs(isc)%cdirname = ALLSPACE_DIR
            myconfs(isc)%cdirname = TRIM(tempdirname)

            filename = ALLSPACE_FILE
            filename = TRIM(ADJUSTL(myconfs(isc)%cfilename))
            myconfs(isc)%cfilename = ALLSPACE_FILE
            myconfs(isc)%cfilename = TRIM(filename)

            filename = ALLSPACE_FILE
            filename = TRIM(MDHEADDIR)//'/'//TRIM(myconfs(isc)%cdirname)//'/'//TRIM(myconfs(isc)%cfilename)
            OPEN(NEWUNIT=CONFSCAN,FILE=TRIM(filename),IOSTAT=ios,STATUS='OLD',ACTION='READ')
            IF(ios .NE. 0) THEN
                WRITE(*,'(A,A)') "ERROR) CANNOT OPEN CONF FILE : ", TRIM(filename)
                STOP
            END IF

            WRITE(*,'(A,A)') " READCONFS) READING FROM CONF FILE : ",TRIM(filename)
            WRITE(LOGUNIT,'(A,A)') " READCONFS) READING FROM CONF FILE : ",TRIM(filename)
            
            ! READ THE NUMBER OF CONFIGURATIONS IN THIS SCAN
            READ(CONFSCAN,*,IOSTAT=ios) 
            READ(CONFSCAN,*,IOSTAT=ios) myconfs(isc)%numconfs
            ALLOCATE(myconfs(isc)%confname(myconfs(isc)%numconfs),myconfs(isc)%cdist(myconfs(isc)%numconfs))
            ALLOCATE(myconfs(isc)%conffilename(myconfs(isc)%numconfs),myconfs(isc)%topolname(myconfs(isc)%numconfs))

            DO i = 1, myconfs(isc)%numconfs
                ! READ CONFIGURATION NAME AND THE CORRESPONDING DISTANCE
                READ(CONFSCAN,*,IOSTAT=ios) myconfs(isc)%confname(i), myconfs(isc)%cdist(i),&
                                            myconfs(isc)%conffilename(i), myconfs(isc)%topolname(i)
                IF(ios .NE. 0) THEN
                    WRITE(*,'(A,A)') "ERROR) WHILE READING CONFS FILE : ", TRIM(filename)
                    STOP
                END IF
            END DO
            CLOSE(CONFSCAN)
        END DO
        CLOSE(CONFUNIT)
    END SUBROUTINE READCONFS

    SUBROUTINE SORTCONFS(myconfs)
        TYPE(conftype), INTENT(INOUT) :: myconfs(:)
        INTEGER :: isc, i, j, imintemp
        REAL(KIND=kind_real) :: rmintemp
        CHARACTER(LEN=len_filename) :: cmintemp

        DO isc = 1, NUMSCANS
            IF( .NOT. ALLOCATED(myconfs(isc)%confname)) STOP "INTERNAL ERROR"
            IF( .NOT. ALLOCATED(myconfs(isc)%cdist)) STOP "INTERNAL ERROR"
            IF( .NOT. ALLOCATED(myconfs(isc)%conffilename)) STOP "INTERNAL ERROR"
            IF( .NOT. ALLOCATED(myconfs(isc)%topolname)) STOP "INTERNAL ERROR"

            DO i = 1, myconfs(isc)%numconfs-1
                rmintemp = myconfs(isc)%cdist(i) ; imintemp = i
                DO j = i+1, myconfs(isc)%numconfs
                    IF(myconfs(isc)%cdist(j) .LT. rmintemp) THEN
                        rmintemp = myconfs(isc)%cdist(j)
                        imintemp = j
                    END IF
                END DO
                IF(imintemp .NE. i) THEN
                    ! SWAP DISTANCE FIRST AND IN THE SAME ORDER SWAP WEIGHTS AND CONF NAMES
                    rmintemp = myconfs(isc)%cdist(imintemp)
                    myconfs(isc)%cdist(imintemp) = myconfs(isc)%cdist(i)
                    myconfs(isc)%cdist(i) = rmintemp
                    ! RAW TRANSFER OF STRINGS 
                    cmintemp = ALLSPACE_FILE ; cmintemp = myconfs(isc)%confname(imintemp)
                    myconfs(isc)%confname(imintemp) = myconfs(isc)%confname(i)
                    myconfs(isc)%confname(i) = ALLSPACE_FILE ; myconfs(isc)%confname(i) = ADJUSTL(cmintemp)
                    ! RAW TRANSFER OF STRINGS 
                    cmintemp = ALLSPACE_FILE ; cmintemp = myconfs(isc)%conffilename(imintemp)
                    myconfs(isc)%conffilename(imintemp) = myconfs(isc)%conffilename(i)
                    myconfs(isc)%conffilename(i) = ALLSPACE_FILE ; myconfs(isc)%conffilename(i) = ADJUSTL(cmintemp)
                    ! RAW TRANSFER OF STRINGS 
                    cmintemp = ALLSPACE_FILE ; cmintemp = myconfs(isc)%topolname(imintemp)
                    myconfs(isc)%topolname(imintemp) = myconfs(isc)%topolname(i)
                    myconfs(isc)%topolname(i) = ALLSPACE_FILE ; myconfs(isc)%topolname(i) = ADJUSTL(cmintemp)
                END IF
            END DO
        END DO
    END SUBROUTINE SORTCONFS

    SUBROUTINE CHECK_QM_CONF_DIST(myqmdata,myconfs)
        TYPE(qmdatatype), INTENT(INOUT) :: myqmdata(:)
        TYPE(conftype), INTENT(INOUT)   :: myconfs(:)
        INTEGER :: isc, i
        REAL(KIND=kind_real) :: rmintemp

        DO isc = 1, NUMSCANS
            IF( .NOT. ALLOCATED(myconfs(isc)%cdist)) STOP "INTERNAL ERROR"
            IF( .NOT. ALLOCATED(myconfs(isc)%confname)) STOP "INTERNAL ERROR"
            IF( .NOT. ALLOCATED(qmdata(isc)%qdist)) STOP "INTERNAL ERROR"
            IF( .NOT. ALLOCATED(qmdata(isc)%qenergy)) STOP "INTERNAL ERROR"

            IF(myconfs(isc)%numconfs .NE. qmdata(isc)%numdata) THEN
                WRITE(*,'(A,I0)') "CHECK-QMCONF) ERROR : INCOMPATIBLE QM DATA AND CONFIGURATIONS IN SCAN : ", isc
                WRITE(*,'(A,I0)') "CHECK-QMCONF) ERROR : NUM QM DATA : ", qmdata(isc)%numdata
                WRITE(*,'(A,I0)') "CHECK-QMCONF) ERROR : NUM CONFS   : ", myconfs(isc)%numconfs
                WRITE(LOGUNIT,'(A,I0)') "CHECK-QMCONF) ERROR : INCOMPATIBLE QM DATA AND CONFIGURATIONS IN SCAN : ", isc
                WRITE(LOGUNIT,'(A,I0)') "CHECK-QMCONF) ERROR : NUM QM DATA : ", qmdata(isc)%numdata
                WRITE(LOGUNIT,'(A,I0)') "CHECK-QMCONF) ERROR : NUM CONFS   : ", myconfs(isc)%numconfs
                STOP
            END IF

            rmintemp = 0.D0
            DO i = 1, myconfs(isc)%numconfs
                IF(ABS(myconfs(isc)%cdist(i)-qmdata(isc)%qdist(i)) .GT. 1000.D0*EPSILON(rmintemp)) THEN
                     WRITE(*,'(A,I0)') "CHECK-QMCONF) ERROR : INCOMPATIBLE QM DATA AND CONFIGURATIONS IN SCAN : ", isc
                     WRITE(*,'(A,F8.4,A,F8.4)') "CHECK-QMCONF) ERROR: THE CONF DISTANCE : ", myconfs(isc)%cdist(i), &
                                                " DOESN'T MATCH QM DISTANCE : ", qmdata(isc)%qdist(i) 
                     STOP
                END IF
            END DO
        END DO
    END SUBROUTINE CHECK_QM_CONF_DIST

    SUBROUTINE RUN_MD(parnum,myitppar,mddata,cumerr,mychoice)
        INTEGER, INTENT(IN)               :: parnum
        TYPE(grotype), INTENT(INOUT)      :: myitppar(:)
        TYPE(qmdatatype), INTENT(INOUT)   :: mddata(:)
        REAL(KIND=kind_real), INTENT(OUT) :: cumerr
        INTEGER, INTENT(IN)               :: mychoice
        REAL(KIND=kind_real), ALLOCATABLE :: scancumerr(:)

        !         
        CHARACTER(LEN=len_dirname)              :: temppardir 
        CHARACTER(LEN=len_dirname)              :: tempplotsdir
        CHARACTER(LEN=len_filename)             :: tempmdpfile
        CHARACTER(LEN=len_filename)             :: tempgtpfilename
        TYPE :: scandirtype
            CHARACTER(LEN=len_dirname)              :: dirname
            CHARACTER(LEN=len_dirname), ALLOCATABLE :: tempconfdir(:)
        END TYPE scandirtype
        TYPE(scandirtype), ALLOCATABLE              :: tempscandir(:)

        ! GROMACS EXECUTION RELATED
        CHARACTER(LEN=len_syscommand) :: mycommand
        INTEGER :: EXSTAT, CSTAT, ICHK
        INTEGER :: ENGPLOT, ERRPLOT, WERRPLOT, SCANPLOT
        CHARACTER(LEN=len_syscommand) :: CMSG
        REAL(KIND=kind_real) :: allenergy(ENERGYPARTS), einf
        
        ! BOOK KEEPING VARIABLES
        INTEGER :: isc, i, j
        CHARACTER(LEN=len_dirname) :: tempdir

        ALLOCATE(scancumerr(NUMSCANS)) 

        ! CREATE PAR DIRECTORY
        temppardir = ALLSPACE_DIR
        IF(mychoice == 0) THEN
            WRITE(temppardir,'(A,A,I0,A)') TRIM(RUNDIR), '/eq-', parnum, '/'
        ELSEIF(mychoice == 1) THEN
            WRITE(temppardir,'(A,A,I0,A)') TRIM(RUNDIR), '/ps-', parnum, '/'
        END IF
        CALL CREATEDIR(temppardir) 
        WRITE(LOGUNIT,'(A,A)') "RUNMD) CREATED PARAMETER DIRECTORY : ", TRIM(temppardir)

        ! LOCATE MDP FILE
        tempmdpfile = ALLSPACE_FILE
        tempmdpfile = TRIM(MDHEADDIR)//'/'//TRIM(mdpfilename)

        ! COPY GUESSPAR.GTP FILE INTO EQ-0 DIRECTORY
        IF(parnum .EQ. 0) THEN
            WRITE(LOGUNIT,'(A)') "RUNMD) COPYING GUESS GTP FILE.. "
            tempgtpfilename = ALLSPACE_FILE
            tempgtpfilename = TRIM(MDHEADDIR)//'/'//TRIM(gtpfilename)
            tempdir = ALLSPACE_DIR
            tempdir = TRIM(temppardir)
            CALL FCOPY(tempgtpfilename,temppardir)
        END IF

        ! CREATE PAR.GTP FILE IN PAR DIRECTORY
        WRITE(LOGUNIT,'(A)') "RUNMD) CREATING GTP FILE FOR THIS PARAMETER SET.. "
        tempgtpfilename = ALLSPACE_FILE
        tempgtpfilename = TRIM(temppardir)//TRIM(INCLUDE_GTP_FILE)
        CALL WRITEITP(tempgtpfilename,myitppar)

        ! CREATE PLOTS DIRECTORY
        WRITE(LOGUNIT,*) " CREATING PLOTS DIRECTORY "
        tempdir = ALLSPACE_FILE
        tempdir = TRIM(temppardir)//'plots/'
        CALL CREATEDIR(tempdir)

        tempfile = ALLSPACE_FILE
        tempfile = TRIM(temppardir)//'plots/scan-rmse.dat'

        OPEN(NEWUNIT=SCANPLOT,FILE=tempfile)
        WRITE(SCANPLOT,'(A)') "# SCAN NUMBER      RMSE      SCAN NAME"

        ALLOCATE(tempscandir(NUMSCANS))
        DO isc = 1, NUMSCANS

            ! CREATE SCAN DIRECTORY
            tempscandir(isc)%dirname = ALLSPACE_DIR
            tempscandir(isc)%dirname = TRIM(temppardir)//TRIM(myconfs(isc)%cscanname)//'/'
            
            WRITE(LOGUNIT,'(A,A)') "RUNMD) CREATING SCAN DIRECTORY : ", TRIM(tempscandir(isc)%dirname)
            CALL CREATEDIR(tempscandir(isc)%dirname)

            ! CREATE CONF DIRECTORIES
            WRITE(LOGUNIT,'(A)') "RUNMD) CREATING CONF SUBDIRECTORIES.. "
            IF(ALLOCATED(tempscandir(isc)%tempconfdir)) DEALLOCATE(tempscandir(isc)%tempconfdir)
            ALLOCATE(tempscandir(isc)%tempconfdir(myconfs(isc)%numconfs))
            DO i = 1, myconfs(isc)%numconfs
                tempscandir(isc)%tempconfdir(i) = ALLSPACE_FILE
                tempscandir(isc)%tempconfdir(i) = TRIM(tempscandir(isc)%dirname)//TRIM(myconfs(isc)%confname(i))//'/'
                CALL CREATEDIR(tempscandir(isc)%tempconfdir(i))
            END DO

            ! COPY TOPOL/PAR FILES AND MDCONFS TO CONF DIRECTORIES
            WRITE(LOGUNIT,'(A)') "RUNMD) COPYING GRO/PDB FILES TO CONF DIRECTORIES.. "
            DO i = 1, myconfs(isc)%numconfs
                ! WRITE(LOGUNIT,'(A)') "RUNMD) COPYING TOP MDP FILES.. "
                
                tempfile = ALLSPACE_FILE
                tempfile = TRIM(MDHEADDIR)//'/'//TRIM(myconfs(isc)%cdirname)//'/'//&
                           TRIM(myconfs(isc)%topolname(i))
                CALL FCOPY(tempfile,tempscandir(isc)%tempconfdir(i))
                CALL FCOPY(tempgtpfilename,tempscandir(isc)%tempconfdir(i))

                
                tempfile = ALLSPACE_FILE
                tempfile = TRIM(MDHEADDIR)//'/'//TRIM(myconfs(isc)%cdirname)//'/'//&
                           TRIM(myconfs(isc)%conffilename(i))
                CALL FCOPY(tempfile,tempscandir(isc)%tempconfdir(i))
            END DO

            WRITE(LOGUNIT,'(A)') "RUNMD) RUNNING GROMACS.. "
            DO i = 1, myconfs(isc)%numconfs
                ! RUN GROMPP
                ! WRITE(*,'(I0,A)') i, " RUNNING GROMPP "
                mycommand = TRIM(GROMACSEXE)//' grompp '//&
                            ' -f '//TRIM(tempmdpfile)//&
                            ' -p '//TRIM(tempscandir(isc)%tempconfdir(i))//TRIM(myconfs(isc)%topolname(i))//&
                            ' -c '//TRIM(tempscandir(isc)%tempconfdir(i))//TRIM(myconfs(isc)%conffilename(i))//&
                            ' -o '//TRIM(tempscandir(isc)%tempconfdir(i))//TRIM(myconfs(isc)%confname(i))//'.tpr'//&
                            ' -po '//TRIM(tempscandir(isc)%tempconfdir(i))//'mdout.mdp'//&
                            ' -pp '//TRIM(tempscandir(isc)%tempconfdir(i))//'full-topol.top'//&
                            ' 1>'//TRIM(tempscandir(isc)%tempconfdir(i))//'grompp-stdout'//&
                            ' 2>'//TRIM(tempscandir(isc)%tempconfdir(i))//'grompp-stderr'

                CALL EXECUTE_COMMAND_LINE(TRIM(mycommand),EXITSTAT=EXSTAT,CMDSTAT=CSTAT,CMDMSG=CMSG)
                IF(CSTAT .NE. 0) THEN
                    WRITE(*,*) "ERROR) COULD NOT RUN GROMPP.."
                    WRITE(*,*) "-->"
                    WRITE(*,*) "ERROR) COMMAND  : ", TRIM(mycommand)
                    WRITE(*,*) "-->"
                    WRITE(*,*) "ERROR) EXIT MSG : ", TRIM(CMSG)
                    STOP
                END IF
                IF(EXSTAT .NE. 0) THEN
                    WRITE(*,*) "ERROR) GROMPP COMMAND GIVES AN ERROR.."
                    WRITE(*,*) "-->"
                    WRITE(*,*) "ERROR) COMMAND  : ", TRIM(mycommand)
                    STOP
                END IF

                    ! RUN MDRUN
                    ! WRITE(*,'(I0,A)') i, " RUNNING MDRUN "
                mycommand = TRIM(GROMACSEXE)//' mdrun '//&
                            ' -s '//TRIM(tempscandir(isc)%tempconfdir(i))//TRIM(myconfs(isc)%confname(i))//'.tpr'//&
                            ' -deffnm '//TRIM(tempscandir(isc)%tempconfdir(i))//TRIM(myconfs(isc)%confname(i))//&
                            ' -rerun '//TRIM(tempscandir(isc)%tempconfdir(i))//TRIM(myconfs(isc)%conffilename(i))//&
                            ' -ntomp '//' 1 '//&
                            ' 1>'//TRIM(tempscandir(isc)%tempconfdir(i))//'mdrun-stdout'//&
                            ' 2>'//TRIM(tempscandir(isc)%tempconfdir(i))//'mdrun-stderr'
                CALL EXECUTE_COMMAND_LINE(TRIM(mycommand),EXITSTAT=EXSTAT,CMDSTAT=CSTAT,CMDMSG=CMSG)
                IF(CSTAT .NE. 0) THEN
                    WRITE(*,*) "ERROR) COULD NOT RUN MDRUN.."
                    WRITE(*,*) "-->"
                    WRITE(*,*) "ERROR) COMMAND  : ", TRIM(mycommand)
                    WRITE(*,*) "-->"
                    WRITE(*,*) "ERROR) EXIT MSG : ", TRIM(CMSG)
                    STOP
                END IF
                IF(EXSTAT .NE. 0) THEN
                    WRITE(*,*) "ERROR) MDRUN GIVES AND ERROR.."
                    WRITE(*,*) "-->"
                    WRITE(*,*) "ERROR) COMMAND  : ", TRIM(mycommand)
                    STOP
                END IF

                    ! RUN ENERGY
                    ! WRITE(*,'(I0,A)') i, " RUNNING GMX ENERGY "
                mycommand = TRIM(GROMACSEXE)//' energy '//&
                            ' -f '//TRIM(tempscandir(isc)%tempconfdir(i))//TRIM(myconfs(isc)%confname(i))//'.edr'//&
                            ' -o '//TRIM(tempscandir(isc)%tempconfdir(i))//'all-energies.xvg'//&
                            ' 1>'//TRIM(tempscandir(isc)%tempconfdir(i))//'energy-stdout'//&
                            ' 2>'//TRIM(tempscandir(isc)%tempconfdir(i))//'energy-stderr'//&
                            ' <<EOF '//mynewline//&
                            ' 1 2 3 4 5 6 7 8 '//&
                            'EOF'//mynewline
                        
                CALL EXECUTE_COMMAND_LINE(TRIM(mycommand),EXITSTAT=EXSTAT,CMDSTAT=CSTAT,CMDMSG=CMSG)
                IF(CSTAT .NE. 0) THEN
                    WRITE(*,*) "ERROR) COULD NOT RUN GMX ENERGY.."
                    WRITE(*,*) "-->"
                    WRITE(*,*) "ERROR) COMMAND  : ", TRIM(mycommand)
                    WRITE(*,*) "-->"
                    WRITE(*,*) "ERROR) EXIT MSG : ", TRIM(CMSG)
                    STOP
                END IF
                IF(EXSTAT .NE. 0) THEN
                    WRITE(*,*) "ERROR) GMX ENERGY GIVES AND ERROR.."
                    WRITE(*,*) "-->"
                    WRITE(*,*) "ERROR) COMMAND  : ", TRIM(mycommand)
                    STOP
                END IF
            END DO

            ! EXTRACT TOTAL ENERGIES
            WRITE(LOGUNIT,'(A)') "RUNMD) EXTRACTING TOTAL ENERGIES "
            DO i = 1, myconfs(isc)%numconfs
                tempdir = ALLSPACE_FILE
                tempdir = TRIM(tempscandir(isc)%tempconfdir(i))//'all-energies.xvg'
                CALL READENERGY(tempdir,allenergy)
                      ! SUBTRACTING BOND, ANGLE, DIHEDRAL FROM PE
                mddata(isc)%qenergy(i) = allenergy(IPOTENERGY) - allenergy(IBONDENERGY) - allenergy(IANGLEENERGY) &
                                         - allenergy(IDIHEDENERGY)
            END DO
            
            einf = mddata(isc)%qenergy(mddata(isc)%numdata)
            mddata(isc)%qenergy(:) = mddata(isc)%qenergy(:) - einf

            WRITE(LOGUNIT,'(A)') "RUNMD) CREATING PLOTS FILES "
            
            tempfile = ALLSPACE_FILE
            tempfile = TRIM(temppardir)//'plots/'//TRIM(myconfs(isc)%cscanname)//'-compare-energy-scans.dat'
            
            OPEN(NEWUNIT=ENGPLOT,FILE=TRIM(tempfile))
            WRITE(ENGPLOT,'(A,I0,A)') "# DIST    QMENERGY    MDENERGY  (ps-",parnum,")"
            DO i = 1, qmdata(isc)%numdata
                WRITE(ENGPLOT,'(F8.4,2X,F16.8,2X,F16.8)') qmdata(isc)%qdist(i), qmdata(isc)%qenergy(i), mddata(isc)%qenergy(i)
            END DO

            tempfile = ALLSPACE_FILE
            tempfile = TRIM(temppardir)//'plots/'//TRIM(myconfs(isc)%cscanname)//'-sqerr.dat'
            
            OPEN(NEWUNIT=ERRPLOT,FILE=tempfile)
            WRITE(ERRPLOT,'(A,I0,A)') "#   DIST        ERROR          ERROR**2  (ps-",parnum,")"
            DO i = 1, qmdata(isc)%numdata
                WRITE(ERRPLOT,'(F8.4,2X,F16.8,2X,F24.8)') qmdata(isc)%qdist(i), mddata(isc)%qenergy(i)-qmdata(isc)%qenergy(i), &
                                                          (mddata(isc)%qenergy(i)-qmdata(isc)%qenergy(i))**2 
            END DO
        
            cumerr = 0.D0 ; scancumerr(isc) = 0.D0
            
            tempfile = ALLSPACE_FILE
            tempfile = TRIM(temppardir)//'plots/'//TRIM(myconfs(isc)%cscanname)//'-wsqerr.dat'
            
            OPEN(NEWUNIT=WERRPLOT,FILE=tempfile)
            WRITE(WERRPLOT,'(A,I0,A)') "# DIST        WEIGHTED ERROR    WEIGHTED ERROR**2  (ps-",parnum,")"
            DO i = 1, qmdata(isc)%numdata
                WRITE(WERRPLOT,'(F8.4,2X,F16.8,2X,F24.8)') qmdata(isc)%qdist(i), &
                                qmdata(isc)%qweight(i)*(mddata(isc)%qenergy(i)-qmdata(isc)%qenergy(i)), &
                                qmdata(isc)%qweight(i)*(mddata(isc)%qenergy(i)-qmdata(isc)%qenergy(i))**2 
                cumerr = cumerr + qmdata(isc)%qweight(i)*(mddata(isc)%qenergy(i)-qmdata(isc)%qenergy(i))**2
            END DO
            cumerr = SQRT(cumerr/SUM(qmdata(isc)%qweight(:)))
            scancumerr(isc) = cumerr

            WRITE(SCANPLOT,'(I2,2X,F16.8,2X,A)') isc, scancumerr(isc), TRIM(myconfs(isc)%cscanname)

            CLOSE(ENGPLOT) ; CLOSE(ERRPLOT) ; CLOSE(WERRPLOT)
        END DO
        CLOSE(SCANPLOT)
        cumerr = SUM(scancumerr(:))/REAL(NUMSCANS)
    END SUBROUTINE RUN_MD

    SUBROUTINE INIT_MC_LOGS(mychoice)
        INTEGER, INTENT(IN) :: mychoice
        CHARACTER(LEN=len_filename) :: tempfile
        CHARACTER(LEN=len_record)   :: fheader, emsg

      ! 0 --> EQUIL ;; 1 --> PRODUCTION
      IF(mychoice == 0) THEN
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!********************************* ALL PARAMETER LOGS **************************************!!    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        tempfile = TRIM(RUNDIR)//'/equil-all-error.dat'
        OPEN(NEWUNIT=MCALLENER,FILE=TRIM(tempfile),IOSTAT=ios)
        emsg = "ERROR) CANNOT OPEN MC LOG FILE : "//TRIM(tempfile) ; CALL CHECKIOS(ios,emsg)
        
        WRITE(MCALLENER,'(A,A)') "# MCITER  ", "    WERROR    "

        tempfile = ALLSPACE_FILE ; tempfile = TRIM(RUNDIR)//'/equil-all-par.dat'
        OPEN(NEWUNIT=MCALLPAR,FILE=TRIM(tempfile),IOSTAT=ios)
        emsg = ALLSPACE_REC ; emsg = "ERROR) CANNOT OPEN MC LOG FILE : "//TRIM(tempfile) ; CALL CHECKIOS(ios,emsg)
      
        fheader = "# MCITER "
        DO i = 1, NUMLJPAR
            fheader = TRIM(fheader)//REPEAT(' ',4)//firstpar(i)%parname(1:4)//"-- SIG"
            fheader = TRIM(fheader)//REPEAT(' ',4)//firstpar(i)%parname(1:4)//"-- EPS"
        END DO
        WRITE(MCALLPAR,'(A)') TRIM(fheader)
        !*********************************************************************************************!

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!******************************* ACCEPTED PARAMETER LOGS ***********************************!!    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        tempfile = ALLSPACE_FILE ; tempfile = TRIM(RUNDIR)//'/equil-acc-error.dat'
        OPEN(NEWUNIT=MCACCENER,FILE=TRIM(tempfile),IOSTAT=ios)
        emsg = ALLSPACE_REC ; emsg = "ERROR) CANNOT OPEN MC LOG FILE : "//TRIM(tempfile) ; CALL CHECKIOS(ios,emsg)

        WRITE(MCACCENER,'(A,A)') "# ACCITER ", "    WERROR    "
        

        tempfile = ALLSPACE_FILE ; tempfile = TRIM(RUNDIR)//'/equil-acc-par.dat'
        OPEN(NEWUNIT=MCACCPAR,FILE=TRIM(tempfile),IOSTAT=ios)
        emsg = ALLSPACE_REC ; emsg = "ERROR) CANNOT OPEN MC LOG FILE : "//TRIM(tempfile) ; CALL CHECKIOS(ios,emsg)
        
        fheader = REPEAT(' ',len_record)
        fheader = "# ACCITER  TOTITER"
        DO i = 1, NUMLJPAR
            fheader = TRIM(fheader)//REPEAT(' ',4)//firstpar(i)%parname(1:4)//"-- SIG"
            fheader = TRIM(fheader)//REPEAT(' ',4)//firstpar(i)%parname(1:4)//"-- EPS"
        END DO
        WRITE(MCACCPAR,'(A)') TRIM(fheader)
        !*********************************************************************************************!

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!******************************* REJECTED PARAMETER LOGS ***********************************!!        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        tempfile = ALLSPACE_FILE ; tempfile = TRIM(RUNDIR)//'/equil-rej-error.dat'
        OPEN(NEWUNIT=MCREJENER,FILE=TRIM(tempfile),IOSTAT=ios)
        emsg = ALLSPACE_REC ; emsg = "ERROR) CANNOT OPEN MC LOG FILE : "//TRIM(tempfile) ; CALL CHECKIOS(ios,emsg)

        WRITE(MCREJENER,'(A,A)') "# REJITER ", "    WERROR    "
        
        tempfile = ALLSPACE_FILE ; tempfile = TRIM(RUNDIR)//'/equil-rej-par.dat'
        OPEN(NEWUNIT=MCREJPAR,FILE=TRIM(tempfile),IOSTAT=ios)
        emsg = ALLSPACE_REC ; emsg = "ERROR) CANNOT OPEN MC LOG FILE : "//TRIM(tempfile) ; CALL CHECKIOS(ios,emsg)

        fheader = REPEAT(' ',len_record)
        fheader = "# REJITER  TOTITER"
        DO i = 1, NUMLJPAR
            fheader = TRIM(fheader)//REPEAT(' ',4)//firstpar(i)%parname(1:4)//"-- SIG"
            fheader = TRIM(fheader)//REPEAT(' ',4)//firstpar(i)%parname(1:4)//"-- EPS"
        END DO
        WRITE(MCREJPAR,'(A)') TRIM(fheader)
        !*********************************************************************************************!
        WRITE(*,'(A)') "MCLOGS) OPENED EQUIL LOG FILES.. "
        WRITE(LOGUNIT,'(A)') "MCLOGS) OPENED EQUIL LOG FILES.. "
      ELSE
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!********************************* ALL PARAMETER LOGS **************************************!!    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        tempfile = TRIM(RUNDIR)//'/all-error.dat'
        OPEN(NEWUNIT=MCALLENER,FILE=TRIM(tempfile),IOSTAT=ios)
        emsg = "ERROR) CANNOT OPEN MC LOG FILE : "//TRIM(tempfile) ; CALL CHECKIOS(ios,emsg)
        
        WRITE(MCALLENER,'(A,A)') "# MCITER  ", "    WERROR    "

        tempfile = ALLSPACE_FILE ; tempfile = TRIM(RUNDIR)//'/all-par.dat'
        OPEN(NEWUNIT=MCALLPAR,FILE=TRIM(tempfile),IOSTAT=ios)
        emsg = ALLSPACE_REC ; emsg = "ERROR) CANNOT OPEN MC LOG FILE : "//TRIM(tempfile) ; CALL CHECKIOS(ios,emsg)
      
        fheader = "# MCITER  "
        DO i = 1, NUMLJPAR
            fheader = TRIM(fheader)//REPEAT(' ',4)//firstpar(i)%parname(1:4)//"-- SIG"
            fheader = TRIM(fheader)//REPEAT(' ',4)//firstpar(i)%parname(1:4)//"-- EPS"
        END DO
        WRITE(MCALLPAR,'(A)') TRIM(fheader)
        !*********************************************************************************************!

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!******************************* ACCEPTED PARAMETER LOGS ***********************************!!    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        tempfile = ALLSPACE_FILE ; tempfile = TRIM(RUNDIR)//'/acc-error.dat'
        OPEN(NEWUNIT=MCACCENER,FILE=TRIM(tempfile),IOSTAT=ios)
        emsg = ALLSPACE_REC ; emsg = "ERROR) CANNOT OPEN MC LOG FILE : "//TRIM(tempfile) ; CALL CHECKIOS(ios,emsg)

        WRITE(MCACCENER,'(A,A)') "# ACCITER ", "    WERROR    "
        

        tempfile = ALLSPACE_FILE ; tempfile = TRIM(RUNDIR)//'/acc-par.dat'
        OPEN(NEWUNIT=MCACCPAR,FILE=TRIM(tempfile),IOSTAT=ios)
        emsg = ALLSPACE_REC ; emsg = "ERROR) CANNOT OPEN MC LOG FILE : "//TRIM(tempfile) ; CALL CHECKIOS(ios,emsg)
        
        fheader = REPEAT(' ',len_record)
        fheader = "# ACCITER  TOTITER"
        DO i = 1, NUMLJPAR
            fheader = TRIM(fheader)//REPEAT(' ',4)//firstpar(i)%parname(1:4)//"-- SIG"
            fheader = TRIM(fheader)//REPEAT(' ',4)//firstpar(i)%parname(1:4)//"-- EPS"
        END DO
        WRITE(MCACCPAR,'(A)') TRIM(fheader)
        !*********************************************************************************************!

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!******************************* REJECTED PARAMETER LOGS ***********************************!!        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        tempfile = ALLSPACE_FILE ; tempfile = TRIM(RUNDIR)//'/rej-error.dat'
        OPEN(NEWUNIT=MCREJENER,FILE=TRIM(tempfile),IOSTAT=ios)
        emsg = ALLSPACE_REC ; emsg = "ERROR) CANNOT OPEN MC LOG FILE : "//TRIM(tempfile) ; CALL CHECKIOS(ios,emsg)

        WRITE(MCREJENER,'(A,A)') "# REJITER ", "    WERROR    "
        
        tempfile = ALLSPACE_FILE ; tempfile = TRIM(RUNDIR)//'/rej-par.dat'
        OPEN(NEWUNIT=MCREJPAR,FILE=TRIM(tempfile),IOSTAT=ios)
        emsg = ALLSPACE_REC ; emsg = "ERROR) CANNOT OPEN MC LOG FILE : "//TRIM(tempfile) ; CALL CHECKIOS(ios,emsg)

        fheader = REPEAT(' ',len_record)
        fheader = "# REJITER  TOTITER"
        DO i = 1, NUMLJPAR
            fheader = TRIM(fheader)//REPEAT(' ',4)//firstpar(i)%parname(1:4)//"-- SIG"
            fheader = TRIM(fheader)//REPEAT(' ',4)//firstpar(i)%parname(1:4)//"-- EPS"
        END DO
        WRITE(MCREJPAR,'(A)') TRIM(fheader)
        !*********************************************************************************************!
        WRITE(*,'(A)') "MCLOGS) OPENED PROD LOG FILES.. "
        WRITE(LOGUNIT,'(A)') "MCLOGS) OPENED PROD LOG FILES.. "
      END IF
    END SUBROUTINE INIT_MC_LOGS

    SUBROUTINE MCMOVE(currljpar,nextljpar,nextitppar)
        TYPE(ljpartype), INTENT(IN)  :: currljpar(:) 
        TYPE(ljpartype), INTENT(OUT) :: nextljpar(:)
        TYPE(grotype), INTENT(OUT)   :: nextitppar(:)
        INTEGER :: i, j

        DO j = 1, NUMLJPAR
            i = 0
            DO 
                i = i + 1
                CALL RANDOM_NUMBER(myrandnum)
                myrandnum = 2.D0*myrandnum - 1.D0
                nextpar(j)%par(:) = currpar(j)%par(:) + currpar(j)%sens(:)*myrandnum(:)
                IF((nextpar(j)%par(1) .GT. nextpar(j)%lwrbnd(1) .AND. nextpar(j)%par(1) .LT. nextpar(j)%uprbnd(1)) &
                .AND. (nextpar(j)%par(2) .GT. nextpar(j)%lwrbnd(2) .AND. nextpar(j)%par(2) .LT. nextpar(j)%uprbnd(2))) THEN
                    EXIT
                END IF
                IF(i .EQ. max_move_tries/2 ) THEN
                    WRITE(*,'(A,A)')       "MCMOVE) WARNING!! TOO MANY TRIES TO GENERATE PAR IN BOUNDS: ",&
                                           TRIM(nextpar(j)%parname)
                    WRITE(LOGUNIT,'(A,A)') "MCMOVE) WARNING!! TOO MANY TRIES TO GENERATE PAR IN BOUNDS: ",&
                                           TRIM(nextpar(j)%parname)
                END IF
                IF(i .EQ. max_move_tries) THEN
                    WRITE(*,'(A,A)')       "MCMOVE) ERROR: TOO MANY TRIES TO GENERATE PAR IN BOUNDS: ",&
                                           TRIM(nextpar(j)%parname)
                    WRITE(LOGUNIT,'(A,A)') "MCMOVE) ERROR: TOO MANY TRIES TO GENERATE PAR IN BOUNDS: ",&
                                           TRIM(nextpar(j)%parname)
                    STOP
                END IF
            END DO
            nextpar(j)%sens(:) = currpar(j)%sens(:)
            curritppar(j)%gropar = currpar(j)
            nextitppar(j)%gropar = nextpar(j)

            WRITE(*,'(A,A,2(2X,F10.6))')       "MCMOVE) --> PAR ", TRIM(nextitppar(j)%gropar%parname),&
                                               nextitppar(j)%gropar%par(:)
            WRITE(LOGUNIT,'(A,A,2(2X,F10.6))') "MCMOVE) --> PAR ", TRIM(nextitppar(j)%gropar%parname),&
                                               nextitppar(j)%gropar%par(:)
        END DO 
    END SUBROUTINE MCMOVE

    SUBROUTINE GETBETA(myiter)
        INTEGER, INTENT(IN) :: myiter
        REAL(KIND=kind_real) :: tempbetaprev, tempbetanext

        IF(myiter .LE. mcbeta%biter(mcbeta%currpnt)) THEN
            mcbeta%bcurr = mcbeta%bcurr + mcbeta%brate
        ELSE
            mcbeta%currpnt = mcbeta%currpnt + 1
            mcbeta%brate = (mcbeta%beta(mcbeta%currpnt) - mcbeta%beta(mcbeta%currpnt-1))/&
                    REAL(mcbeta%biter(mcbeta%currpnt) - mcbeta%biter(mcbeta%currpnt-1),KIND=kind_real)
            mcbeta%bcurr = mcbeta%bcurr + mcbeta%brate
            WRITE(*,'(A)')       "GETBETA) BETA UPDATED TO THE NEXT ANNEALING REGIME"
            WRITE(LOGUNIT,'(A)') "GETBETA) BETA UPDATED TO THE NEXT ANNEALING REGIME"
        END IF
    END SUBROUTINE GETBETA

    SUBROUTINE SORTBETA()
        INTEGER :: i, j, minpos, miniter
        REAL(KIND=kind_real) :: minbeta

        DO i = 1 , mcbeta%numpoints-1
            miniter = mcbeta%biter(i)
            minpos = i
            DO j = i+1, mcbeta%numpoints
                IF(miniter .GT. mcbeta%biter(j)) THEN
                    miniter = mcbeta%biter(j)
                    minpos = j
                END IF
            END DO
            IF(minpos .NE. i) THEN
                mcbeta%biter(minpos) = mcbeta%biter(i)
                mcbeta%biter(i) = miniter
                minbeta = mcbeta%beta(minpos)
                mcbeta%beta(minpos) = mcbeta%beta(i)
                mcbeta%beta(i) = minbeta
            END IF
        END DO
        mcbeta%beta(0) = mcbeta%binit ; mcbeta%beta(mcbeta%numpoints+1) = mcbeta%beta(mcbeta%numpoints)
        IF(MAXVAL(mcbeta%biter(1:mcbeta%numpoints)) .GT. MCPROD) THEN
            WRITE(*,'(A,I0,A,I0)') "SORTBETA) ERROR: ANNEALING POINT: ", MAXVAL(mcbeta%biter(1:mcbeta%numpoints)), &
                           " IS GREATER THAN PRODUCITON RUN ITERATIONS: ", MCPROD
            STOP
        END IF
        mcbeta%biter(0) = 0 ; mcbeta%biter(mcbeta%numpoints+1) = MCPROD
        mcbeta%currpnt = 1  ; mcbeta%bcurr = mcbeta%binit
        mcbeta%brate = (mcbeta%beta(1) - mcbeta%beta(0))/REAL(mcbeta%biter(1) - mcbeta%biter(0),KIND=kind_real)
    END SUBROUTINE SORTBETA

    SUBROUTINE INITSENS()
        curr_acc_iter = 0 
        curr_acc_ratio = 0.D0 ;  prev_acc_ratio = 0.D0 ; preverr_acc_ratio = 0.D0
        err_acc_ratio = 0.D0 ; interr_acc_ratio = 0.D0 ; diferr_acc_ratio = 0.D0 
        corr_factor = 0.D0
    END SUBROUTINE INITSENS

    SUBROUTINE INIT_ACC_LOGS()
        CHARACTER(LEN=len_filename) :: tempacclogfile, emsg
        CHARACTER(LEN=len_record)   :: fheader
        INTEGER :: tempios, i

        ! COPY THE INITIAL GUESSPAR.GTP FILE INTO SENS 
        tempacclogfile = TRIM(MDHEADDIR)//'/'//TRIM(gtpfilename)
        CALL FCOPY(tempacclogfile,SENSDIR)

        tempacclogfile = ALLSPACE_FILE
        tempacclogfile = TRIM(RUNDIR)//'acceptance-ratio.dat'
        
        OPEN(NEWUNIT=UNIT_ACCLOG,FILE=TRIM(tempacclogfile),IOSTAT=tempios)
        emsg = "ERROR) CANNOT OPEN ACCEPTANCE RATIO LOG FILE : "//TRIM(tempacclogfile) 
        CALL CHECKIOS(tempios,emsg)

        WRITE(UNIT_ACCLOG,'(A,T2,A,T14,A,T29,A,T44,A,T59,A)') "#", " ITER ", " ACC-RAT ", &
                                                    " DEL-ACC ", " INT-ERR ", " DIF-ERR "

        tempacclogfile = ALLSPACE_FILE
        tempacclogfile = TRIM(RUNDIR)//'sens.dat'
        
        OPEN(NEWUNIT=UNIT_SENSLOG,FILE=TRIM(tempacclogfile),IOSTAT=tempios)
        emsg = "ERROR) CANNOT OPEN SENS LOG FILE : "//TRIM(tempacclogfile) 
        CALL CHECKIOS(tempios,emsg)

        fheader = "# MCITER "
        DO i = 1, NUMLJPAR
            fheader = TRIM(fheader)//REPEAT(' ',4)//firstpar(i)%parname(1:4)//"-- SIG"
            fheader = TRIM(fheader)//REPEAT(' ',4)//firstpar(i)%parname(1:4)//"-- EPS"
        END DO
        WRITE(UNIT_SENSLOG,'(A)') TRIM(fheader)
    END SUBROUTINE INIT_ACC_LOGS

    SUBROUTINE UPDATESENS()
        INTEGER :: i
        corr_factor = kp*err_acc_ratio + ki*interr_acc_ratio + kd*diferr_acc_ratio
        ! AUTO ADJUST THE corr_factor MAGNITUDE TO LESS THAN 1
        DO
            IF(ABS(corr_factor) .LT. 1.D0) EXIT
            corr_factor = corr_factor*0.1D0
            WRITE(*,'(A)')       "UPDATESENS) CORRECTION FACTOR EXCEEDED 1.0 SCALING BACK.."
            WRITE(LOGUNIT,'(A)') "UPDATESENS) CORRECTION FACTOR EXCEEDED 1.0 SCALING BACK.."
        END DO

        ! CURRENT ACC-RAT IS LESS THAN SET VALUE ==> REDUCE SENS(:)
        IF(corr_factor .LT. 0.D0) THEN
            ! RIGHT NOW BOTH KIND (+/-) OF CORRECTIONS ARE APPLIED IN SIMILAR MANNER 
            DO i = 1, NUMLJPAR
                currpar(i)%sens(1) = currpar(i)%sens(1) + sens_scale*firstpar(i)%sens(1)*corr_factor
                IF(currpar(i)%sens(1) .LT. 0) THEN
                    currpar(i)%sens(1) = MINSENS_SIG
                END IF
                currpar(i)%sens(2) = currpar(i)%sens(2) + sens_scale*firstpar(i)%sens(2)*corr_factor
                IF(currpar(i)%sens(2) .LT. 0) THEN
                    currpar(i)%sens(2) = MINSENS_EPS
                END IF
            END DO
        ELSE
            DO i = 1, NUMLJPAR
                currpar(i)%sens(1) = currpar(i)%sens(1) + sens_scale*firstpar(i)%sens(1)*corr_factor
                IF(currpar(i)%sens(1) .LT. 0) THEN
                    currpar(i)%sens(1) = MINSENS_SIG
                END IF
                currpar(i)%sens(2) = currpar(i)%sens(2) + sens_scale*firstpar(i)%sens(2)*corr_factor
                IF(currpar(i)%sens(2) .LT. 0) THEN
                    currpar(i)%sens(2) = MINSENS_EPS
                END IF
            END DO
        END IF
        DO i = 1, NUMLJPAR
            WRITE(*,'(A,A,A,2(2X,F10.6))')       "UPDATESENS)", " --> SENS ", TRIM(currpar(i)%parname),&
                                               currpar(i)%sens(:)
            WRITE(LOGUNIT,'(A,A,A,2(2X,F10.6))') "UPDATESENS)", " --> SENS ", TRIM(currpar(i)%parname),&
                                               currpar(i)%sens(:)
            curritppar(i)%gropar%sens(:) = currpar(i)%sens(:)
        END DO
    END SUBROUTINE UPDATESENS

    SUBROUTINE CHECKIOS(myios,emsg)
        INTEGER, INTENT(IN) :: myios
        CHARACTER(LEN=len_record), INTENT(IN) :: emsg

        IF(ios .NE. 0) THEN
            WRITE(*,'(A)') TRIM(emsg)
            STOP
        END IF
    END SUBROUTINE CHECKIOS
END PROGRAM main
