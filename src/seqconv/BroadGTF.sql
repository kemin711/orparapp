/** load GTF from Broad into a relational database 

	Broad have very simple formats
**/

-- there is no gene table for Broad

drop table if exists transcript_raw;
create table transcript_raw (
	transcriptid varchar(36) primary key,
	geneid varchar(36),
	key gid (geneid)
);
-- not attribute for Broad transcripts
load data local infile '/home/kzhou/work/analysis/externalGenomes/Coprinus_cinereus/transcript.tab'
into table transcript_raw
ignore 1 lines
;

drop table if exists transcript;
create table transcript (
	transcriptid integer primary key,
	geneid integer
);

insert into transcript
select substring_index(transcriptid, '_', -1) as transcriptid,
	substring_index(geneid, '_', -1) as geneid
from transcript_raw;



create table supercontig (
	id integer primary key,
	length integer,
	start_contig integer,
	stop_contig integer
);
load data local infile 'supercontigs.csv'
into table supercontig fields terminated by ','
ignore 1 lines;

drop table if exists contig;
create table contig (
	id INTEGER primary key,
	length integer,
	supercontig integer references supercontig
);
load data local infile 'contigs.csv.txt'
into table contig
fields terminated by ','
ignore 1 lines;

/** do we realy need the exon number?
*/
drop table if exists exon_raw;
create table exon_raw (
	genomicid integer not null references contig, 
	strand char(1), 
	start integer, 
	end integer, 
	transcriptid varchar(36)
);
load data local infile 'exon.tab' into table exon_raw
ignore 1 lines;


drop table if exists exon;
create table exon (
	id serial primary key,
	genomicid integer not null references contig, 
	strand char(1), 
	start integer, 
	end integer, 
	transcriptid integer references transcript, 
	unique key uk (genomicid, transcriptid, start),
	key loc (genomicid, strand, start, end)
);


insert into exon
(genomicid, strand, start, end, transcriptid)
select genomicid, strand, start, end,
	substring_index(transcriptid, '_', -1) as transcriptid
from exon_raw
where strand = '+'
order by genomicid, strand, transcriptid,
	start, end
;

insert into exon
(genomicid, strand, start, end, transcriptid)
select genomicid, strand, start, end,
	substring_index(transcriptid, '_', -1) as transcriptid
from exon_raw
where strand = '-'
order by genomicid, strand, transcriptid,
	start DESC, end DESC
;

select *
from exon
order by genomicid, id;

/** CDS **/
drop table if exists CDS_raw;
create table CDS_raw (
	genomicid integer,
	strand char(1),
	phase integer,
	start integer,
	end integer,
	transcriptid varchar(36)
);

load data local infile 'CDS.tab' into table CDS_raw
ignore 1 lines;


drop table if exists CDS;
create table CDS (
	id serial primary key,
	genomicid integer not null references contig, 
	strand char(1), 
	phase integer,
	start integer, 
	end integer, 
	transcriptid integer references transcript, 
	unique key uk (genomicid, transcriptid, start),
	key loc (genomicid, strand, start, end)
);

insert into CDS
(genomicid, strand, phase, start, end, transcriptid)
select genomicid, strand, phase, start, end,
	substring_index(transcriptid, '_', -1) as transcriptid
from CDS_raw
where strand = '+'
order by genomicid, strand, transcriptid,
	start, end
;

insert into CDS
(genomicid, strand, phase, start, end, transcriptid)
select genomicid, strand, phase, start, end,
	substring_index(transcriptid, '_', -1) as transcriptid
from CDS_raw
where strand = '-'
order by genomicid, strand, transcriptid,
	start DESC, end DESC
;

select *
from CDS
order by genomicid, id;

/** start and stop is associated with CDS information **/
drop table if exists start_codon_raw;
create table start_codon_raw (
	genomicid integer, 
	strand char(1), 
	start integer, end integer, 
	transcriptid varchar(36)
);
load data local infile 'start_codon.tab'
into table start_codon_raw ignore 1 lines;

create table start_codon (
	genomicid integer,
	strand char(1),
	start integer, end integer,
	transcriptid integer references transcript,
	primary key (transcriptid, start)
);

insert into start_codon 
select genomicid, strand, start, end, 
	substring_index(transcriptid, '_', -1) as transcriptid
from start_codon_raw
;

/** stop codon **/
drop table if exists stop_codon_raw;
create table stop_codon_raw (
	genomicid integer, 
	strand char(1), 
	start integer, end integer, 
	transcriptid varchar(36)
);
load data local infile 'stop_codon.tab'
into table stop_codon_raw ignore 1 lines;

create table stop_codon (
	genomicid integer,
	strand char(1),
	start integer, end integer,
	transcriptid integer references transcript,
	primary key (transcriptid, start)
);

insert into stop_codon 
select genomicid, strand, start, end, 
	substring_index(transcriptid, '_', -1) as transcriptid
from stop_codon_raw
;

create table contigseq (
	id integer primary key,
	seq longtext
);


load data local infile 'contigs.tab'
into table contigseq;

/** I guess the proteinid is the same as the transcriptid
*/

drop table if exists proteinseq;
create table protein (
	id integer primary key,
	title varchar(250),
	seq longtext
);

load data local infile '/home/kzhou/work/analysis/externalGenomes/Coprinus_cinereus/proteins.fasta.tab'
into table protein;

select p.id, t.*
from proteinseq p join transcript t
	on p.id=t.transcriptid
;

-- the same number of transcripts and protein


