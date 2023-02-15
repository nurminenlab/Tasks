t1 = tcp('172.27.85.52', 51007);
t2 = tcp('172.27.85.52', 51008);

tcp_write(t1,'a');
tcp_write(t2,num2str(0.5));
pause(1);
tcp_write(t1,'b');
tcp_write(t2,num2str(0.8));
pause(1);

tcp_write(t1,'c');
tcp_write(t2,num2str(0.8));
pause(1);

tcp_write(t1,'b');
tcp_write(t2,num2str(0.6));
pause(1);

tcp_write(t1,'c');
tcp_write(t2,num2str(0.5));
pause(1);

tcp_write(t1,'a');
tcp_write(t2,num2str(0.5));
pause(1);

tcp_write(t1,'b');
tcp_write(t2,num2str(0.8));
tcp_write(t1,'c');
tcp_write(t2,num2str(0.8));
tcp_write(t1,'b');
tcp_write(t2,num2str(0.6));
tcp_write(t1,'c');
tcp_write(t2,num2str(0.5));
tcp_write(t1,'a');
tcp_write(t2,num2str(0.5));
tcp_write(t1,'b');
tcp_write(t2,num2str(0.8));
tcp_write(t1,'c');
tcp_write(t2,num2str(0.8));
tcp_write(t1,'b');
tcp_write(t2,num2str(0.6));
tcp_write(t1,'c');
tcp_write(t2,num2str(0.5));
tcp_write(t1,'a');
tcp_write(t2,num2str(0.5));
tcp_write(t1,'b');
tcp_write(t2,num2str(0.8));
tcp_write(t1,'c');
tcp_write(t2,num2str(0.8));
tcp_write(t1,'b');
tcp_write(t2,num2str(0.6));
tcp_write(t1,'c');
tcp_write(t2,num2str(0.5));
tcp_write(t1,'a');
tcp_write(t2,num2str(0.5));
tcp_write(t1,'b');
tcp_write(t2,num2str(0.8));
tcp_write(t1,'c');
tcp_write(t2,num2str(0.8));
tcp_write(t1,'b');
tcp_write(t2,num2str(0.6));
tcp_write(t1,'c');
tcp_write(t2,num2str(0.5));
tcp_write(t1,'a');
tcp_write(t2,num2str(0.5));
tcp_write(t1,'b');
tcp_write(t2,num2str(0.8));
tcp_write(t1,'c');
tcp_write(t2,num2str(0.8));

# to END
tcp_write(t1,'x');
tcp_write(t2,num2str(0.0));


close t1;
close t2;


