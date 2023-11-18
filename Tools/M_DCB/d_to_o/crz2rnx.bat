@echo off

rem  batch program : CRZ2RNX.bat
rem        created by HATANAKA, Y.  19/Dec/1996
rem        e-mail: hata@gsi-mc.go.jp
rem               (hata@gsi.go.jp after 2001.01.01)
rem  RINEX file compression
rem  *** wildcard can be used ***
rem
rem  compressed CompactRINEX(*.??e)   --- CompactRINEX(*.??d) --- RINEX OBS file (*.??o)
rem  CompactRINEX(*.??d)              --- RINEX OBS file (*.??o)
rem  compressed RINEX NAV file(*.??x) --- RINEX NAV file(*.??n)
rem  compressed RINEX MET file(*.??w) --- RINEX MET file(*.??m)

for %%f in (%1) do call crz2rnx1 %%f

