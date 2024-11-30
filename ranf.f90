     Real(Kind=8)  function  ranf(iq)

       implicit none
       integer iq
       integer IP,IR
       parameter (IP = 48828125, IR = 2147483647)
       
       iq=iq* IP
       !       print *,'iq = ',iq
       if(iq) 10,20,20
10     iq=(iq+IR)+1
20     ranf = dble(iq)/2.0D0**31
       
     end function ranf

