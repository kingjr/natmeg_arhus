
ttl                 = [];
ttl.fields.side     = 1:2;
ttl.fields.target   = 1:6;
ttl.names           = param2ttlnames(cfg.ttl.fields);
ttl.dec             = name2ttl(cfg.ttl.names);
for trigger = 1:12
    disp(ttl2name(trigger,ttl));
end