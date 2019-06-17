from simplesocket import SimpleServer
from datetime import datetime
import sys
import os
import learner


class LearningServer(object):
	def __init__(self):
		assert len(sys.argv) == 3
		self.server = SimpleServer.server_from_string(sys.argv[1])
		runDir = sys.argv[2]
		self.carboxLearner = learner.carboxLearner.CarboxLearner(runDir, "Regression")
		print("Time:{}\t Started:{}".format(datetime.now(), self.server))
		sys.stdout.flush()


	def run(self):
		def clients_handel(client, server):
			cmd = client.receive()
			print("""{} >>> {}""".format(client, cmd))
			assert cmd == "[Join]"
			msg = "{} joined. #Clients : {}".format(client, len(server.clients))
			client.send(msg)
			while True:
				try:
					cmd = client.receive()
					print("""{} >>> {}""".format(client, cmd))
					if cmd != "[Exit]":
						cmd = cmd.split()
						if cmd[0] == "[Time]":
							msg = 'Time is {}'.format(datetime.now().time())
							client.send(msg)
						elif cmd[0] == "[Predict]":
							assert len(cmd) == 2
							msg = str(self.carboxLearner.predict(cmd[1]))
							client.send(msg)
						elif cmd[0] == "[AddData]":
							assert len(cmd) == 2
							msg = str(self.carboxLearner.addData(cmd[1]))
							print(msg)
							client.send(msg)
						else:
							msg = "Error : unknown command."
							client.send(msg)
					else:
						server.clients.remove(client)
						msg = "{} left. #Clients : {}".format(client, len(server.clients))
						client.send(msg)
						client.close()
						break
				except:
					pass

		self.server.run(clients_handel)




if __name__ == '__main__':
	mysever = LearningServer()
	mysever.run()

