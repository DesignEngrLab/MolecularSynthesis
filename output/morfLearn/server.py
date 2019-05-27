from simplesocket import SimpleServer
from datetime import datetime
import sys
import os
import learner


class LearningServer(object):
    def __init__(self):
        # self.server = SimpleServer.server_from_string("" if len(sys.argv) < 2 else sys.argv[1])
        # print("\n\t{}\n\tStarted {}\n".format(datetime.now(), self.server))
        dataDir = os.path.join(os.getcwd(), "..", "tempData")
        self.carboxLearner = learner.carboxLearner.CarboxLearner(dataDir, "Regression")

    def learn(self):
        self.carboxLearner.learn_fix()


    def clients_handel(self, client, server):
        print("{} joined. # clients : {}".format(client, len(server.clients)))
        while True:
            try:
                cmd = client.receive()
                if cmd not in ("[Exit]", "[E]"):
                    print("""{} >>> {}""".format(client, cmd))
                    msg = "Error : unknown command."
                    cmd = cmd.split()
                    if cmd[0] == "time":
                        msg = 'time is {}'.format(datetime.now().time())
                        client.send(msg)
                    elif cmd[0] == "[ALL]":
                        server.broadcast(" ".join(cmd[1:]))
                    else:
                        client.send(msg)
                else:
                    client.close()
                    server.clients.remove(client)
                    print("""{} left. #Clients : {}""".format(client, len(server.clients)))
                    break
            except:
                pass


def main():
    mysever = LearningServer()
    mysever.learn()



if __name__ == '__main__':
    main()
