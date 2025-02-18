from PlasmaProfiles import PlasmaProfiles
import pdb

profStart = PlasmaProfiles(0.57, 0.58, 1.85, 0.7, 0.0, 0.0, 12.188)
profStart.plot_Profiles()

profStable = PlasmaProfiles(0.2, 0.57, 1.85, 0.7, 2.5e6, 9e6, 12.207)
profStable.plot_Profiles()

profPreTQ = PlasmaProfiles(0.35, 0.41, 1.85, 0.7, 2.5e6, 9.87e6, 12.2)
profPreTQ.plot_Profiles()


#profIpSpike = PlasmaProfiles(0.57, 0.58, 0.7, 0.3, 0.0, 4.78e6, 12.104) # 10% Ip spike
profIpSpike = PlasmaProfiles(0.57, 0.58, 1.85, 0.7, 0.0, 5.29e6, 12.085)
profIpSpike.plot_Profiles()


profEnd = PlasmaProfiles(0.57, 0.58, 1.85, 0.7, 0.0, 0.0, 12.188)
profEnd.plot_Profiles()