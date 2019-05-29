program main

use test, only: test_value
use io, only: stdout

use calendar

implicit none

integer :: date(7), ierr
character(len=128) :: string
double precision :: jd

string = '2001-01-01'
call parse_date(string,date,'YYYY-MM-DD',ierr)
call test_value(date(1),2001,'parse_date: YYYY-MM-DD: YYYY')
call test_value(date(2),1,'parse_date: YYYY-MM-DD: MM')
call test_value(date(3),1,'parse_date: YYYY-MM-DD: DD')
call test_value(date(4),0,'parse_date: YYYY-MM-DD: HH')
call test_value(date(5),0,'parse_date: YYYY-MM-DD: MM')
call test_value(date(6),0,'parse_date: YYYY-MM-DD: SS')
call test_value(date(7),0,'parse_date: YYYY-MM-DD: MMM')

string = '1992 02 29'
call parse_date(string,date,'YYYY MM DD',ierr)
call test_value(date(1),1992,'parse_date: YYYY MM DD: YYYY')
call test_value(date(2),2,'parse_date: YYYY MM DD: MM')
call test_value(date(3),29,'parse_date: YYYY MM DD: DD')
call test_value(date(4),0,'parse_date: YYYY MM DD: HH')
call test_value(date(5),0,'parse_date: YYYY MM DD: MM')
call test_value(date(6),0,'parse_date: YYYY MM DD: SS')
call test_value(date(7),0,'parse_date: YYYY MM DD: MMM')

string = '19891225'
call parse_date(string,date,'YYYYMMDD',ierr)
call test_value(date(1),1989,'parse_date: YYYYMMDD: YYYY')
call test_value(date(2),12,'parse_date: YYYYMMDD: MM')
call test_value(date(3),25,'parse_date: YYYYMMDD: DD')
call test_value(date(4),0,'parse_date: YYYYMMDD: HH')
call test_value(date(5),0,'parse_date: YYYYMMDD: MM')
call test_value(date(6),0,'parse_date: YYYYMMDD: SS')
call test_value(date(7),0,'parse_date: YYYYMMDD: MMM')

string = '2010-02-27T11:12:13'
call parse_date(string,date,'YYYY-MM-DDTHH:MM:SS',ierr)
call test_value(date(1),2010,'parse_date: YYYY-MM-DDTHH:MM:SS: YYYY')
call test_value(date(2),2,'parse_date: YYYY-MM-DDTHH:MM:SS: MM')
call test_value(date(3),27,'parse_date: YYYY-MM-DDTHH:MM:SS: DD')
call test_value(date(4),11,'parse_date: YYYY-MM-DDTHH:MM:SS: HH')
call test_value(date(5),12,'parse_date: YYYY-MM-DDTHH:MM:SS: MM')
call test_value(date(6),13,'parse_date: YYYY-MM-DDTHH:MM:SS: SS')
call test_value(date(7),0,'parse_date: YYYY-MM-DDTHH:MM:SS: MMM')

string = '2010 02 27 11 12 13'
call parse_date(string,date,'YYYY MM DD HH MM SS',ierr)
call test_value(date(1),2010,'parse_date: YYYY MM DD HH MM SS: YYYY')
call test_value(date(2),2,'parse_date: YYYY MM DD HH MM SS: MM')
call test_value(date(3),27,'parse_date: YYYY MM DD HH MM SS: DD')
call test_value(date(4),11,'parse_date: YYYY MM DD HH MM SS: HH')
call test_value(date(5),12,'parse_date: YYYY MM DD HH MM SS: MM')
call test_value(date(6),13,'parse_date: YYYY MM DD HH MM SS: SS')
call test_value(date(7),0,'parse_date: YYYY MM DD HH MM SS: MMM')

string = '20100227111213'
call parse_date(string,date,'YYYYMMDDHHMMSS',ierr)
call test_value(date(1),2010,'parse_date: YYYYMMDDHHMMSS: YYYY')
call test_value(date(2),2,'parse_date: YYYYMMDDHHMMSS: MM')
call test_value(date(3),27,'parse_date: YYYYMMDDHHMMSS: DD')
call test_value(date(4),11,'parse_date: YYYYMMDDHHMMSS: HH')
call test_value(date(5),12,'parse_date: YYYYMMDDHHMMSS: MM')
call test_value(date(6),13,'parse_date: YYYYMMDDHHMMSS: SS')
call test_value(date(7),0,'parse_date: YYYYMMDDHHMMSS: MMM')

string = '2011-03-11T14:15:16.567'
call parse_date(string,date,'YYYY-MM-DDTHH:MM:SS.MMM',ierr)
call test_value(date(1),2011,'parse_date: YYYY-MM-DDTHH:MM:SS.MMM: YYYY')
call test_value(date(2),3,'parse_date: YYYY-MM-DDTHH:MM:SS.MMM: MM')
call test_value(date(3),11,'parse_date: YYYY-MM-DDTHH:MM:SS.MMM: DD')
call test_value(date(4),14,'parse_date: YYYY-MM-DDTHH:MM:SS.MMM: HH')
call test_value(date(5),15,'parse_date: YYYY-MM-DDTHH:MM:SS.MMM: MM')
call test_value(date(6),16,'parse_date: YYYY-MM-DDTHH:MM:SS.MMM: SS')
call test_value(date(7),567,'parse_date: YYYY-MM-DDTHH:MM:SS.MMM: MMM')

string = '2011 03 11 14 15 16 567'
call parse_date(string,date,'YYYY MM DD HH MM SS MMM',ierr)
call test_value(date(1),2011,'parse_date: YYYY MM DD HH MM SS MMM: YYYY')
call test_value(date(2),3,'parse_date: YYYY MM DD HH MM SS MMM: MM')
call test_value(date(3),11,'parse_date: YYYY MM DD HH MM SS MMM: DD')
call test_value(date(4),14,'parse_date: YYYY MM DD HH MM SS MMM: HH')
call test_value(date(5),15,'parse_date: YYYY MM DD HH MM SS MMM: MM')
call test_value(date(6),16,'parse_date: YYYY MM DD HH MM SS MMM: SS')
call test_value(date(7),567,'parse_date: YYYY MM DD HH MM SS MMM: MMM')

string = '20110311141516567'
call parse_date(string,date,'YYYYMMDDHHMMSSMMM',ierr)
call test_value(date(1),2011,'parse_date: YYYYMMDDHHMMSSMMM: YYYY')
call test_value(date(2),3,'parse_date: YYYYMMDDHHMMSSMMM: MM')
call test_value(date(3),11,'parse_date: YYYYMMDDHHMMSSMMM: DD')
call test_value(date(4),14,'parse_date: YYYYMMDDHHMMSSMMM: HH')
call test_value(date(5),15,'parse_date: YYYYMMDDHHMMSSMMM: MM')
call test_value(date(6),16,'parse_date: YYYYMMDDHHMMSSMMM: SS')
call test_value(date(7),567,'parse_date: YYYYMMDDHHMMSSMMM: MMM')


call date2jd(date,jd)
call test_value(jd,2455632.0939417477d0,'date2jd')
date = 0
call jd2date(jd,date)
call test_value(date(1),2011,'jd2date: year')
call test_value(date(2),3,'jd2date: month')
call test_value(date(3),11,'jd2date: day')
call test_value(date(4),14,'jd2date: hour')
call test_value(date(5),15,'jd2date: min')
call test_value(date(6),16,'jd2date: sec')
! call test_value(date(7),567,'jd2date: ms')


write(stdout,*) 'calendar_module unit test passed'
end
