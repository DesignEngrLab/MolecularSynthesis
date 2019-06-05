from simplesocket import SimpleServer
from datetime import datetime
import sys
import os
import learner


class LearningServer(object):
	def __init__(self):
		self.server = SimpleServer.server_from_string("" if len(sys.argv) < 2 else sys.argv[1])
		self.carboxLearner = learner.carboxLearner.CarboxLearner(os.path.join(os.getcwd(), "computation", "data"), "Regression")
		print("Time:{}\t Started:{}".format(datetime.now(), self.server))
		sys.stdout.flush()


	def run(self):
		def clients_handel(client, server):
			msg = "{} joined. # clients : {}".format(client, len(server.clients))
			client.send(msg)
			while True:
				try:
					cmd = client.receive()
					if cmd != "[Exit]":
						print("""{} >>> {}""".format(client, cmd))
						cmd = cmd.split()
						if cmd[0] == "[Time]":
							msg = 'Time is {}'.format(datetime.now().time())
							client.send(msg)
						elif cmd[0] == "[Predict]":
							msg = self.carboxLearner.predict(cmd[1]).__str__()
							client.send(msg)
						else:
							msg = "Error : unknown command."
							client.send(msg)
					else:
						server.clients.remove(client)
						msg = """{} left. #Clients : {}""".format(client, len(server.clients))
						client.send(msg)
						client.close()
						break
				except:
					pass

		self.server.run(clients_handel)




if __name__ == '__main__':
	mysever = LearningServer()
	mysever.run()

