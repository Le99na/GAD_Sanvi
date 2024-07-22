Hi Sanvi,

here some seeds to start your EEG journey:

- functional script on GAD data --> TF_MC.m
	-phaseAmp_calculator.m
	-pinkFreqScale.m


- on Dropbox I have found Justins Code for the 2PAC study (see paper)
	-pac_2PAC_Rest_sept2020.m
	-pac_2PAC_Task_sept2020.m
	-pac_2PAC-Task_sept2020_tgdb.m
	--> they are quite similar, have not yet figured out the differences


- I have also downloaded a script time-frequency in 2PAC (timeFreq_2PAC_Task_aug2020.m) and about wPLI in 2 PAC (wPLI_2PAC_Rest_feb2021.m) which I am not sure I you need them.


The following steps should be first in order not to loose oversight and get lost in the code
1) look at our GAD data and think about what time within the task would be relevant to look at (where would cognitive control be most prevalent)
2) what is the format of our data and how would you define the conditions (justin has always so many more conditions than we have)
	- right now we have the format: channelsxtimextrials
3) then you can look at Justin's code and try to decide which parts are relevant for us and which are not
	