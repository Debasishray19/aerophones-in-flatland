# Aerophones In Flatland - Test Procedure And Test Result

<b >Running Code</b>.
<br>1. Main File: aeophonesInFlatLand_v2.m
<br>2. To inject impulse as source use the following MATLAB file: <i>impulseResponse.m</i> and comment the following code in  <i>aeophonesInFlatLand_v2.m </i>
```diff
- excitationV = srcAmplitude * sin(2*pi*excitationF*dt*(exeT(:)-1));
```
<br>3. To verify the frequency response of audio data inlude the MATLAB file: <i>audioGenfunc.m</i> and run the following code line
```diff
+ audioGenfunc(Pr_Audio, excitationV);
```


