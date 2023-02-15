t1 = tcp('172.27.85.52', 51009);
disp("connected");
tcp_write(t1,100);
tcp_write(t1,100);
tcp_write(t1,100);
tcp_write(t1,100);
tcp_write(t1,num2str(100));
clear t1;
