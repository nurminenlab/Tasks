import socket
import time
j = "Plotter Ready!"
s1 = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s1.bind(('172.27.85.105', 51006))
s1.listen(1)
conn, addr = s1.accept()
print("connected from ", addr,"to recieve stim_cond")
'''conn.sendall(f"{j}".encode())
#time.sleep(0.001)
j1 = "Plotter Ready to recieve runtime"

s2 = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s2.bind(('localhost', 51006))
s2.listen(1)
conn2, addr2 = s2.accept()
print("connected from ", addr2,"to runtime for INTAN")
conn2.sendall(f"{j1}".encode())'''
#time.sleep(0.001)

while True :
    msg = conn.recv(1024).decode()
 #   runtime = conn2.recv(1024).decode()

    # receive data stream. it won't accept data packet greater than 1024 bytes
    if not msg or msg == 'end':
        # if data is not received break
        break
    print(msg )
#    print(int(runtime))

conn.close()
#conn2.close()