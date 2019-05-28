from simplesocket import SimpleServer
from datetime import datetime
import sys
import os
import learner


class LearningServer(object):
    def __init__(self):
        self.server = SimpleServer.server_from_string("" if len(sys.argv) < 2 else sys.argv[1])
        # print("\n\t{}\n\tStarted {}\n".format(datetime.now(), self.server))
        # self.carboxLearner = learner.carboxLearner.CarboxLearner(os.path.join(os.getcwd(), "computation", "data"))

    # def learn(self):
    #     self.carboxLearner.learn()

    def run(self):
        def clients_handel(client, server):
            print("{} joined. # clients : {}".format(client, len(server.clients)))
            while True:
                try:
                    cmd = client.receive()
                    if cmd != "[E]":
                        print("""{} >>> {}""".format(client, cmd))
                        cmd = cmd.split()
                        if cmd[0] == "time":
                            msg = 'time is {}'.format(datetime.now().time())
                            client.send(msg)
                        else:
                            msg = "Error : unknown command."
                            client.send(msg)
                    else:
                        client.close()
                        server.clients.remove(client)
                        print("""{} left. #Clients : {}""".format(client, len(server.clients)))
                        break
                except:
                    pass

        self.server.run(clients_handel)




def main():
    mysever = LearningServer()
    mysever.run()



if __name__ == '__main__':
    sys.argv[1] = "9993"
    main()
