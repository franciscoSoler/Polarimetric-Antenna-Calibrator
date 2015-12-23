
max_gain = 0;
index = 0;
for rm = 11:20
	for i = 1:32
		pass = false;
		for j = 1:10
			if wal2(j, :) == ph_wal(i, :)
				pass = true;
				break
			end
		end
		if pass
			continue
		end
		wal2(rm, :) = ph_wal(i,:);
		[g, p] = pccdecoder(in(1:amountPulses:length(in(:,1)), :), chirp.chirpreplica(10,:), wal2);
		if g(rm) > 54
			rm 
			index = i
			max_gain = g(rm)
		end
	end
end