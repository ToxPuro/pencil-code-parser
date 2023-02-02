some code here

! and comments here

!$omp parallel shared(f,df) !$parser-command:generate-firstprivate-pragma-for-static-variables
!$omp do
