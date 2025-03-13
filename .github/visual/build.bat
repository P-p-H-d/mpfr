REM It shall be run from ROOTSRC
REM Overwrite ports of VCPKG ro reuse its build infrastructure
REM with latest version of MPFR
xcopy /y .github\visual\tst.patch C:\vcpkg\ports\mpfr
xcopy /y .github\visual\portfile.cmake C:\vcpkg\ports\mpfr
REM Run install of mpfr by VCPKG
cd C:\vcpkg
REM Disable failing test with VC++
set MPFR_TESTS_SKIP_CHECK_NULL=1
vcpkg install mpfr:x64-windows --head
set "MyReturnCode=%ERRORLEVEL%"
REM Display outputs of build
type C:\vcpkg\buildtrees\mpfr\config-x64-windows-dbg-out.log
type C:\vcpkg\buildtrees\mpfr\config-x64-windows-dbg-err.log
type C:\vcpkg\buildtrees\mpfr\build-x64-windows-dbg-out.log
type C:\vcpkg\buildtrees\mpfr\build-x64-windows-dbg-err.log
exit /B %MyReturnCode%
