program MOPAC_win 
  USE IFQWIN
  use dfport
  use molkst_C, only: jobnam, gui, ijulian, program_name, verson, site_no, line
  use chanel_C, only : iw0
  implicit none
  INTEGER :: i4  
  integer i, j, k, loop
  character :: ch
    gui = .false.
    iw0 = 0
    ijulian = 123
    site_no = 88
    program_name = "Standalone MOPAC "
    OPEN(0, file='user', title='MOPAC 2016') 
    i4 = setbkcolor(15)
    i4 = settextcolor(0)
    call clearscreen($GCLEARSCREEN)
    i = iargc()
    if (i == 0) then
      write(0,*)
      write(0,*)"  MOPAC in WINDOWS: To run MOPAC 2016 within WINDOWS, see:"      
      write(0,*)
      write(0,*)"  http://OpenMOPAC.net/manual/Running_MOPAC.html"
      write(0,*)
      write(0,*)
      write(0,*)
      i = 0
      stop
    end if
    call getdatestamp(line, verson)  
    call password
    call run_mopac()  ! the user has supplied data sets as arguments
    i = SETEXITQQ(QWIN$EXITNOPERSIST)
  end program MOPAC_win
  
  
!+++++++++++++++++++++++++++++++++++++++++++++++++
!  InitialSettings Do not edit, unless necessary
!
LOGICAL FUNCTION InitialSettings()

USE DFLIB
USE DFNLS
implicit none
LOGICAL :: l4, result
integer(4) cp, i
  type (qwinfo) :: qwi
external :: instructions, about, shutdown

call NLSGetLocale(CODEPAGE = cp)   
! This routine is called automatically when the program begins.  It sets
! up the menu structure for the program, and connects "callback" routines
! with each menu item.
  qwi.x = 0
  qwi.y = 0
  qwi.w = 1000
  qwi.h = 1000
  qwi.type = QWIN$SET
  i = setWSizeQQ(QWIN$FRAMEWINDOW, qwi)
   l4 = appendmenuqq(1, $MENUENABLED,   'File'C,           NUL)
   l4 = appendmenuqq(1, $MENUENABLED,   'Exit'C,           WINEXIT)

   l4 = appendmenuqq(2, $MENUENABLED,   'Edit'C,	   NUL)
   l4 = appendmenuqq(2, $MENUENABLED,   'Select Text'C,	   WINSELECTTEXT)
   l4 = appendmenuqq(2, $MENUENABLED,   'Select All'C,	   WINSELECTALL)
   l4 = appendmenuqq(2, $MENUENABLED,   'Copy'C,	   WINCOPY)
   l4 = appendmenuqq(2, $MENUENABLED,   'Paste'C,	   WINPASTE)

   l4 = appendmenuqq(3, $MENUENABLED,   'View'C,           NUL)
   l4 = appendmenuqq(3, $MENUENABLED,   'Size to Fit'C,    WINSIZETOFIT)
   l4 = appendmenuqq(3, $MENUENABLED,   'Full Screen'C,    WINFULLSCREEN)

   l4 = appendmenuqq(4, $MENUENABLED,   'State'C,          NUL)
   l4 = appendmenuqq(4, $MENUENABLED,   'Pause'C,          WINSTATE)

   l4 = appendmenuqq(5, $MENUENABLED,   'Window'C,         NUL)
   l4 = appendmenuqq(5, $MENUENABLED,   'Status Bar'C,     WINSTATUS)

   l4 = appendmenuqq(6, $MENUENABLED,   'Help'C,           NUL)
   l4 = appendmenuqq(6, $MENUENABLED,   'Instructions'C, Instructions)
   l4 = appendmenuqq(6, $MENUENABLED,   'About'C,               About)
   
   l4 = appendmenuqq(7, $MENUENABLED,   'Shutdown'C,         NUL)
   l4 = appendmenuqq(7, $MENUENABLED,   'Shutdown'C,         Shutdown)


InitialSettings = l4

return
  end
  
  subroutine Instructions(checked)
      USE DFLIB
      USE DFNLS
      logical :: checked, lstatus
      integer :: i4
      lstatus = checked
      i4 = messageboxqq('                    Instructions for running MOPAC for WINDOWS in screen modes\r\r &
      & MOPAC for WINDOWS runs MOPAC2016 in screen and command-prompt modes.\r &
      & These instructions apply to the screen mode only.\r\r &
      & "Drag and Drop" mode:\r &
      & Locate the data set to be run, and drag it to the MOPAC for WINDOWS icon.\r\r &
      & The data set will be run, and the results placed in the folder &
      & where the data set came from.\r\r &
      & "Open with" mode:\r &
      & Within Windows Explorer, locate the file to be run\r &
      & Right-click the file, and use "Open with" -> "MOPAC2016.exe"\r\r &
      & If, at that point,  you can''t see "MOPAC2016.exe" then:\r &
      & Use "Choose program" -> "Browse" then navigate to the file "MOPAC2016.exe".\r &
      & Highlight it, then click on "Open".\r &
      & 'C,'MOPAC Instructions'C, MB$OK)
      return
end
subroutine About(checked)
      USE DFLIB
      USE DFNLS
      logical :: checked, lstatus
      integer :: i4
      lstatus = checked
      i4 = messageboxqq('   MOPAC for WINDOWS\r\r &
      & Copyright: Stewart Computational Chemistry\r\r &
      & Version 2016 (2016)\r\r &
      &  'C,'MOPAC Instructions'C, MB$OK)
      return
end
subroutine Shutdown(checked)
      USE DFLIB
      USE DFNLS
      use chanel_C, only: iend, end_fn
      logical :: checked
      open (unit=iend, file=end_fn, status="UNKNOWN", err=1000)
      write(iend,*)"Shut"
      close(iend)
1000  return 
end



