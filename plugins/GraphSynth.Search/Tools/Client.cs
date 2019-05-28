using System;
using System.Net;
using System.Net.Sockets;
using System.Text;


namespace GraphSynth.Search.Tools
{
    public class MessageClient
    {

        private readonly int messagePort;
        private readonly Socket sender;


        public MessageClient(int port)
        {
            messagePort = port;
        }

        public void Connect()
        {
            // Establish the remote endpoint for the socket and connect to the port given
            var hostName = Dns.GetHostName();
            IPHostEntry ipHost = Dns.GetHostEntry(hostName);

            //Here ipHost contains two ip address, only the second one matches python server
            var ipAddr = ipHost.AddressList[1];
            IPEndPoint localEndPoint = new IPEndPoint(ipAddr, messagePort);

            sender = new Socket(ipAddr.AddressFamily, SocketType.Stream, ProtocolType.Tcp);

            // Connect Socket to the remote endpoint using method Connect() 
            try
            {
                Console.WriteLine("Try to connect to server.....");
                sender.Connect(localEndPoint);

                // We print EndPoint information  
                // that we are connected 
                Console.WriteLine("Socket connected to -> {0} ", sender.RemoteEndPoint.ToString());
            }
            catch (Exception e)
            {
                Console.WriteLine(e.ToString());
            }
        }

        public void SendMessage(string msg)
        {
            try
            {
                // Creation of messagge that 
                // we will send to Server 
                byte[] messageSent = Encoding.ASCII.GetBytes(msg);
                int byteSent = sender.Send(messageSent);

                // Data buffer 
                byte[] messageReceived = new byte[1024];

                // We receive the messagge using the method Receive(). 
                // This method returns number of bytes received, 
                // that we'll use to convert them to string
                int byteRecv = sender.Receive(messageReceived);
                Console.WriteLine("Message from Server -> {0}", Encoding.ASCII.GetString(messageReceived, 0, byteRecv));

                // Close Socket using the method Close() 
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


    }


}