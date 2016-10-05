CC 	= gcc
channel7: channel7.c beamformer.c signal.c noise.c IIR_filter.c FIR_filter.c
	$(CC) -g -o $@ $^ -lm -DFLAT

clean:
	$(RM)  channel7

.PHONY: clean
