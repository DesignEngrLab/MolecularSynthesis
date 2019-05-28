using System;
using System.Net;
using System.Net.Sockets;
using System.Text;


namespace GraphSynth.Search.Tools
{
    public class MessageClient
    {

        private readonly int portUsed;

        public MessageClient(int port)
        {
            portUsed = port;
        }

        public void SendMessage(string msg)
        {
            try
            {
                // Establish the remote endpoint  
                // for the socket. This example  
                // uses port 11111 on the local  
                // computer. 
                var hostName = Dns.GetHostName();
                Console.WriteLine(hostName);
                IPHostEntry ipHost = Dns.GetHostEntry(hostName);
                Console.WriteLine(ipHost);
                //Here ipHost contains two ip address, only the second one matches python server.
                var ipAddr = ipHost.AddressList[1];
                IPEndPoint localEndPoint = new IPEndPoint(ipAddr, portUsed);
                Console.WriteLine(localEndPoint);

                // Creation TCP/IP Socket using  
                // Socket Class Costructor 
                Socket sender = new Socket(ipAddr.AddressFamily, SocketType.Stream, ProtocolType.Tcp);
                Console.WriteLine("Try to send message");

                try
                {

                    // Connect Socket to the remote  
                    // endpoint using method Connect() 
                    sender.Connect(localEndPoint);

                    // We print EndPoint information  
                    // that we are connected 
                    Console.WriteLine("Socket connected to -> {0} ",
                                  sender.RemoteEndPoint.ToString());

                    // Creation of messagge that 
                    // we will send to Server 
                    byte[] messageSent = Encoding.ASCII.GetBytes(msg);
                    int byteSent = sender.Send(messageSent);

                    // Data buffer 
                    byte[] messageReceived = new byte[1024];

                    // We receive the messagge using  
                    // the method Receive(). This  
                    // method returns number of bytes 
                    // received, that we'll use to  
                    // convert them to string 
                    int byteRecv = sender.Receive(messageReceived);
                    Console.WriteLine("Message from Server -> {0}",
                          Encoding.ASCII.GetString(messageReceived,
                                                     0, byteRecv));

                    // Close Socket using  
                    // the method Close() 
                    sender.Shutdown(SocketShutdown.Both);
                    sender.Close();
                }

                // Manage of Socket's Exceptions 
                catch (ArgumentNullException ane)
                {
                    Console.WriteLine("ArgumentNullException : {0}", ane.ToString());
                }

                catch (SocketException se)
                {
                    Console.WriteLine("SocketException : {0}", se.ToString());
                }

                catch (Exception e)
                {
                    Console.WriteLine("Unexpected exception : {0}", e.ToString());
                }
            }

            catch (Exception e)
            {
                Console.WriteLine(e.ToString());
            }
        }


    }


}