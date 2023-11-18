@echo off

rem  batch program : CRZ2RNX1.bat 
rem        created by HATANAKA, Y.  19/Dec/1996
rem        e-mail: hata@gsi-mc.go.jp
rem        (hata@gsi.go.jp after 2001.01.01)
rem  RINEX file recovery
rem  uncompress and/or recover RINEX file
rem  *** Only 1 file can be processed ***
rem
rem  compressed CompactRINEX(*.??e)   --- CompactRINEX(*.??d) --- RINEX OBS file (*.??o)
rem  CompactRINEX(*.??d)              --- RINEX OBS file (*.??o)
rem  compressed RINEX NAV file(*.??x) --- RINEX NAV file(*.??n)
rem  compressed RINEX MET file(*.??w) --- RINEX MET file(*.??n)
rem
rem  %1 filename

splname %1 > SPL_NAME.BAT
call SPL_NAME.BAT
del  SPL_NAME.BAT

if %last%==d goto CRX
if %last%==D goto CRX
if %last%==e goto CRZ
if %last%==E goto CRZ
if %last%==y goto ZRX
if %last%==Y goto ZRX
if %last%==x goto NAV
if %last%==X goto NAV
if %last%==w goto MET
if %last%==W goto MET
goto EXIT

:CRX
    crx2rnx %1 -f
    goto EXIT

:CRZ
    copy %base%e %base%z > nul:
    call decompr %base%z %base%d
    crx2rnx %base%d -f
    del %base%z
    del %base%d
    goto EXIT

:ZRX
    copy %base%y %base%z > nul:
    call decompr %base%z %base%o
    del %base%z
    goto EXIT

:NAV
    copy %base%x %base%z > nul:
    call decompr %base%z %base%n
    del  %base%z
    goto EXIT

:MET
    copy %base%w %base%z > nul:
    call decompr %base%z %base%m
    del  %base%z
    goto EXIT

:EXIT
