using System;
using System.Net;
using System.Net.Sockets;
using System.Text;
using System.Threading;


namespace GraphSynth.Search.Tools
{
    public class MessageClient
    {

        private readonly int messagePort;
        private Socket sender;
        private readonly IPEndPoint localEndPoint;

        private static Mutex sendMessageMutex = new Mutex();


        public MessageClient(int port)
        {
            messagePort = port;

            // Establish the remote endpoint for the socket and connect to the port given
            IPHostEntry ipHost = Dns.GetHostEntry(Dns.GetHostName());

            //Here ipHost contains two ip address, only the second one matches python server
            localEndPoint = new IPEndPoint(ipHost.AddressList[1], messagePort);

        }

        public void Connect()
        {
            sender = new Socket(AddressFamily.InterNetwork, SocketType.Stream, ProtocolType.Tcp);
            // Connect Socket to the remote endpoint using method Connect() 
            try
            {
                Console.WriteLine("Try to connect to server.....");
                sender.Connect(localEndPoint);

                // We print EndPoint information  
                // that we are connected
                Console.WriteLine("Socket connected to -> {0} ", sender.RemoteEndPoint.ToString());

                SendMessage("[Join]");
            }
            catch (Exception e)
            {
                Console.WriteLine(e.ToString());
            }
        }

        public void DisConnect()
        {
            SendMessage("[Exit]");
            // Close Socket using the method Close() 
            sender.Shutdown(SocketShutdown.Both);
            sender.Close();
        }

        public string SendMessage(string msg)
        {
            string response = null;
            try
            {
                // Creation of messagge that 
                // we will send to Server 
                byte[] messageSent = Encoding.ASCII.GetBytes(msg);

                // Data buffer 
                byte[] messageReceived = new byte[1024];

                sendMessageMutex.WaitOne();
                sender.Send(messageSent);

                // We receive the messagge using the method Receive(). 
                // This method returns number of bytes received, 
                // that we'll use to convert them to string
                int byteRecv = sender.Receive(messageReceived);
                sendMessageMutex.ReleaseMutex();


                response = Encoding.ASCII.GetString(messageReceived, 0, byteRecv);
                //Console.WriteLine("Message from Server -> {0}", response);
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
            return response;
        }


    }


}