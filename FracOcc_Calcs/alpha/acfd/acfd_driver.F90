subroutine acfd_driver
  use memory, only: memory_save, memory_release
  character(64) :: cmd

  ibase = memory_save()
  call inputl(ncl)
  call gets(1, cmd, 1)
  if (cmd .eq. 'RIRPA') then
    call acfd_rirpa
  elseif (cmd .eq. 'URIRPA') then
    call acfd_urirpa
  elseif (cmd .eq. 'SCEXX') then
    call acfd_scexx
  elseif (cmd .eq. 'USCEXX') then
    call acfd_uscexx
  elseif (cmd .eq. 'KSINV') then
    call ksinv_driver(cmd)
  elseif (cmd .eq. 'UKSINV') then
    call uksinv_driver(cmd)
  end if

  call memory_release(ibase)
  return
end subroutine acfd_driver
