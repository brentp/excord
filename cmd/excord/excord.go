// This program extracts a bedpe of discordant and split reads.
package main

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"log"
	"os"
	"sort"
	"strconv"
	"strings"
	"sync"

	arg "github.com/alexflint/go-arg"
	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"github.com/brentp/bigly"
	"github.com/brentp/bigly/bamat"
	"github.com/brentp/faidx"
	"github.com/brentp/goleft/covmed"
	"github.com/brentp/mmslice/uint16mm"
	"github.com/brentp/windowmean"
	"github.com/brentp/xopen"
)

const minAlignSize = 30

type cliarg struct {
	Prefix             string  `arg:"-p,help:prefix for reference files. if not given they are not written"`
	ExcludeFlag        uint16  `arg:"-F"`
	MinMappingQuality  uint8   `arg:"-Q"`
	DiscordantDistance int     `arg:"-d,help:distance at which mates are considered discordant. if not provided it is calcuated from data"`
	Fasta              string  `arg:"-f,help:path to fasta file. used to check for mismatches"`
	BamPath            string  `arg:"positional,required"`
	Region             string  `arg:"positional"`
	medianReadLength   float64 `arg:"-"`
}

type pairType int

const (
	idiscordantSA = -1
	idiscordant   = 0
)

type bedPE struct {
	c1      string
	s1, e1  int
	strand1 int8

	c2      string
	s2, e2  int
	strand2 int8

	iType int
}

type excord struct {
	readCov            *uint16mm.Slice
	pairCov            *uint16mm.Slice
	altMask            []bool
	discordantDistance int

	chrom string

	f  io.Writer
	ch chan *bedPE
	wg *sync.WaitGroup
}

func cap255(s []uint16) {
	for i, v := range s {
		if v == 255 {
			s[i] = 254
		} else if v > 255 {
			s[i] = 255
		}
	}
}

// modified from http://stackoverflow.com/a/4073700
func round(num uint16, factor uint16) uint16 {
	v := num + factor/2 - (num-1)%factor
	if factor&1 != 0 {
		return v
	}
	return v - 1
}

func _quantize(val uint16) uint16 {
	if val < 7 {
		return val
	}
	if val < 18 {
		if val&1 != 0 {
			return val - 1
		}
		return val
	}
	if val < 48 {
		return round(val, 5)
	}
	if val < 84 {
		return round(val, 9)
	}
	val = round(val, 17)
	if val > 255 {
		return 255
	}
	return val
}

func (ex *excord) quantize() {
	log.Println("quantizing")
	if ex.altMask == nil {
		return
	}
	radius := 7
	ex.pairCov.Flush()
	ex.readCov.Flush()
	cap255(ex.pairCov.A)
	cap255(ex.readCov.A)
	ex.pairCov.Flush()
	ex.readCov.Flush()

	rc, err := uint16mm.Open(nil, len(ex.readCov.A))
	if err != nil {
		panic(err)
	}
	copy(rc.A, ex.readCov.A)
	tmp := windowmean.WindowMeanUint16(rc.A, radius)
	copy(rc.A, tmp)

	pc, err := uint16mm.Open(nil, len(ex.pairCov.A))
	if err != nil {
		panic(err)
	}
	copy(pc.A, ex.pairCov.A)
	tmp = windowmean.WindowMeanUint16(pc.A, radius)
	copy(pc.A, tmp)

	ex.pairCov.Flush()
	ex.readCov.Flush()

	orc, opc := ex.readCov.A, ex.pairCov.A

	for i, m := range ex.altMask {
		// no quantization if we're in an alt region
		_ = m
		/*
			if m {
				continue
			}
		*/
		orc[i] = _quantize(rc.A[i])
		opc[i] = _quantize(pc.A[i])
	}
	log.Println("done quantizing")
}

func newExcord(chromLen int, prefix string, discordantDistance int) *excord {
	e := &excord{altMask: make([]bool, chromLen)}
	rcov, err := uint16mm.Create(prefix+"read.bin", int64(chromLen))
	pcheck(err)
	pcov, err := uint16mm.Create(prefix+"pair.bin", int64(chromLen))
	pcheck(err)
	e.altMask = make([]bool, int64(chromLen))
	e.readCov = rcov
	e.pairCov = pcov
	e.ch = make(chan *bedPE, 5)
	e.wg = &sync.WaitGroup{}
	e.wg.Add(1)
	e.discordantDistance = discordantDistance
	e.f = bufio.NewWriter(os.Stdout)
	go func() {
		for b := range e.ch {
			if _, err := fmt.Fprintf(e.f, "%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\n", b.c1, b.s1, b.e1, b.strand1, b.c2, b.s2, b.e2, b.strand2, b.iType); err != nil {
				panic(err)
			}
			e.updateMask(b)
		}
		e.wg.Done()
	}()

	return e
}

func (ex *excord) updateMask(b *bedPE) {
	if ex.discordantDistance == 0 {
		return
	}
	mask := ex.altMask
	var s, e int
	if b.iType == idiscordant || b.iType == idiscordantSA {
		if b.c1 == ex.chrom {
			if b.strand1 == 1 {
				s, e = b.e1-10, b.e1+ex.discordantDistance
			} else {
				s, e = b.s1-ex.discordantDistance, b.s1+10
			}
			if s < 0 {
				s = 0
			}
			if e > len(mask) {
				e = len(mask)
			}
			for i := s; i < e; i++ {
				mask[i] = true
			}
		}

		if b.c2 == ex.chrom {
			if b.strand2 == 1 {
				s, e = b.e2-10, b.e2+ex.discordantDistance
			} else {
				s, e = b.s2-ex.discordantDistance, b.s2+10
			}
			if s < 0 {
				s = 0
			}
			if e > len(mask) {
				e = len(mask)
			}
			for i := s; i < e; i++ {
				mask[i] = true
			}
		}

		return
	}
	// ##################################
	// splitter
	// ##################################
	// we just want to be able to check the actual split location so we goleft
	// +/- 10 of the observed split.
	// we know that c1:s1-e1 is the left end and c2:s2-e2 is the right
	//  [xxxxxx]-----------------[xxxxxx]
	//       *****             *****
	if b.c1 == ex.chrom {
		s, e := b.e1-10, b.e1+10
		if s < 0 {
			s = 0
		}
		if e > len(mask) {
			e = len(mask)
		}
		for i := s; i < e; i++ {
			mask[i] = true
		}
	}

	if b.c2 == ex.chrom {
		s, e := b.s2-10, b.s2+10
		if s < 0 {
			s = 0
		}
		if e > len(mask) {
			e = len(mask)
		}
		for i := s; i < e; i++ {
			mask[i] = true
		}
	}
}

func (e *excord) Close() error {
	if e.wg != nil {
		e.wg.Wait()
	}
	if b, ok := e.f.(*bufio.Writer); ok {
		b.Flush()
	}
	if e.readCov != nil {
		e.readCov.Close()
		return e.pairCov.Close()
	}
	return nil
}

func (e *excord) WriteAlt(b *bedPE) {
	e.ch <- b
}

func (e *excord) updatePairCoverage(start, end int) {
	A := e.pairCov.A
	for p := end - 1; p >= start; p-- {
		A[p]++
	}
}

func (e *excord) updateReadCoverage(start, end int) {
	A := e.readCov.A
	for p := end - 1; p >= start; p-- {
		A[p]++
	}
}

func (c cliarg) Version() string {
	return "excord 0.2.0"
}

func pcheck(e error) {
	if e != nil {
		panic(e)
	}
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func iabs(a int) int {
	if a < 0 {
		return -a
	}
	return a
}

type interval struct {
	chrom []byte
	start int
	end   int
}

func intStrand(b bool) int8 {
	if b {
		return 1
	}
	return -1
}

// called from writeDiscordant. We've collected the reads and SA tags into slices, now go through all cases and check for problems.
func writeSAs(chrom string, lefts, rights []*bigly.SA, ex *excord, discordantDistance int) []*interval {
	var refs []*interval
	for _, l := range lefts {
		for _, r := range rights {
			cmp := bytes.Compare(l.Chrom, r.Chrom)
			if cmp > 0 || cmp == 0 && r.Pos < l.Pos {
				l, r = r, l
			}
			if discordantSA(l, r, discordantDistance) {
				b := bedPE{string(stripChr(l.Chrom)), l.Pos, l.End(), intStrand(l.Strand), string(stripChr(r.Chrom)), r.Pos, r.End(), intStrand(r.Strand), -1}
				ex.WriteAlt(&b)
			} else if cmp == 0 && bytes.Equal(r.Chrom, []byte(chrom)) {
				refs = append(refs, &interval{r.Chrom, l.Pos, r.End()})

			}
		}
	}
	return refs
}

func writeRefIntervals(refs []*interval, ex *excord) {
	for _, i := range refs {
		ex.updatePairCoverage(i.start, i.end)
	}
}

func discordantSA(a, b *bigly.SA, discordantDistance int) bool {
	var p1, p2 int
	if a.Pos < b.Pos {
		p1, p2 = a.Pos, b.End()
	} else {
		p1, p2 = b.Pos, a.End()
	}
	d := iabs(p1 - p2)
	return d > discordantDistance
}

func getMateEnd(r *sam.Record, opts *cliarg) int {
	if mc, ok := r.Tag([]byte{'M', 'C'}); ok {
		cig, err := sam.ParseCigar(mc[3:])
		if err == nil {

			rl, _ := cig.Lengths()
			return r.MatePos + rl
		}
	}
	return r.MatePos + int(opts.medianReadLength)
}

const bad = sam.Duplicate | sam.Secondary | sam.Supplementary

// output discordant mates to file. The map keeps a history of left mates (lower position) that we've
// seen that have SA tags. We put the read, and it's SA tags into the map.
// Reads without SA tags don't need to go into the map because we can get their start and end from the mate.
// With a pair like:
// [@R1(a)...]------------[@R1(b)....]  [@R2......]
// the left-most splitter (a) of R1 will be a discordant pair; the right-most (b) will be concordant.
// We want to track both pieces of evidence. Once we get to the primay chunk of R1, we put it and it's
// tags into the map.
// Once we get to R2 then we look for its tags, for it's mate in the map.
// We do the outer product of all pairs of R1 and its SA tags and all of R2 and it's SA tags. Some pairs
// may be concordant and some discordant. We return concordant intervals from this function to be increment the
// reference counts.
// Only uses: paired, primary, non-dup, non-supplement reads
func writeDiscordant(r *sam.Record, ex *excord, opts *cliarg, m map[string][]*bigly.SA) (refs []*interval) {
	// we are Paired, not dup or secondary and we are the right-most read on the same chrom.
	if r.Flags&bad != 0 || r.Flags&sam.Paired != sam.Paired {
		return
	}
	if r.Start() == r.MatePos && r.Ref.ID() == r.MateRef.ID() {
		return
	}

	if r.Start() > r.MatePos {
		// we find a discordant pair directly.
		if r.Ref.ID() != r.MateRef.ID() || discordantByDistance(r, opts.DiscordantDistance) {
			start := r.Start()
			mateStart, mateEnd := r.MatePos, getMateEnd(r, opts)
			mateFlag := int8(1)
			if r.Flags&sam.MateReverse == sam.MateReverse {
				mateFlag = -1
			}
			chrom, mateChrom := r.Ref.ID(), r.MateRef.ID()
			// always output left-most mate first.
			if chrom < mateChrom || chrom == mateChrom && start < mateStart {
				b := bedPE{sstripChr(r.Ref.Name()), start, r.End(), r.Strand(), sstripChr(r.MateRef.Name()), mateStart, mateEnd, mateFlag, 0}
				ex.WriteAlt(&b)
			} else {
				b := bedPE{sstripChr(r.MateRef.Name()), mateStart, mateEnd, mateFlag, sstripChr(r.Ref.Name()), start, r.End(), r.Strand(), 0}
				ex.WriteAlt(&b)
			}
		}
		// for each tag in this pair we see if it is
		var asas, bsas []*bigly.SA
		if pt, ok := r.Tag([]byte{'S', 'A'}); ok {
			if bytes.Count([]byte(pt), []byte{';'}) < 3 {
				asas = bigly.AsSAs(r, []byte(pt))
			}
		}
		var ok bool
		bsas, ok = m[r.Name]
		// if they are both empty, we can return.
		if !ok && len(asas) == 0 {
			return
		}
		if ok {
			delete(m, r.Name)
		}
		// otherwise, we create a single SA.
		if len(asas) == 0 {
			asas = []*bigly.SA{bigly.RecordToSA(r)}
		} else if len(bsas) == 0 {
			v, ok := r.Tag([]byte{'M', 'C'})
			if !ok {
				//fmt.Fprintf(os.Stderr, "didn't find mate cigar for %s assuming full match\n", r.Name)
				v = []byte(fmt.Sprintf("MCZ%.0fM", opts.medianReadLength))
			}
			bsas = []*bigly.SA{&bigly.SA{Chrom: []byte(r.MateRef.Name()), Pos: r.MatePos, Cigar: v[3:]}}
		}
		return writeSAs(r.Ref.Name(), asas, bsas, ex, opts.DiscordantDistance)
	}
	if r.MateRef.ID() != r.Ref.ID() {
		return
	}

	// if we are here, we are the left of pair. If we have Tags, add them to the map.
	if t, ok := r.Tag([]byte{'S', 'A'}); ok {
		if bytes.Count([]byte(t), []byte{';'}) <= 3 {
			m[r.Name] = bigly.AsSAs(r, []byte(t))
		}
	}
	return
}

// output splitters. Splitters are ordered by their offset into the read.
// given, cigars of:
// A:20S30M100S
// B:50S30M50S
// C:90S30M30S
// we would order them as they are listed. We would output bedpe intervals
// for A-B, and B-C
func writeSplitter(r *sam.Record, ex *excord) int {
	if r.Flags&sam.Secondary != 0 && r.Flags&sam.Unmapped != 0 {
		return 0
	}
	x, ok := r.Tag([]byte{'S', 'A'})
	if !ok {
		return 0
	}

	var sa []byte = x
	if sa[len(sa)-1] == ';' {
		sa = sa[3 : len(sa)-1]
	} else {
		sa = sa[3:]
	}

	sas := bytes.Split(sa, []byte{';'})
	tags := make([]bigly.SA, len(sas)+1)
	tags[0] = bigly.SA{Chrom: []byte(r.Ref.Name()), Pos: r.Start(), Parsed: r.Cigar}
	tags[0].End()
	// call End() to set Parse internally.
	for i, b := range sas {
		tags[i+1] = bigly.ParseSA(b)
		tags[i+1].End()
	}

	// sort tags by their offset into the query.
	sort.Slice(tags, func(i, j int) bool {
		return bigly.FirstMatch(tags[i].Parsed) < bigly.FirstMatch(tags[j].Parsed)
	})
	// now we have tags sorted by the position in the read that they represent.
	// given splits a, b, c, d we output a -> b, b -> c, c -> d where each ->
	// represents evidence of splitter between those chunks.
	for i, b := range tags {
		if i == 0 {
			continue
		}
		// the last column indicates number of SA fields. Is 0 for discordants above.
		a := tags[i-1]
		cmp := bytes.Compare(a.Chrom, b.Chrom)
		// always output the left-most first.
		if cmp < 0 || cmp == 0 && a.Pos < b.Pos {
			b := bedPE{string(stripChr(a.Chrom)), a.Pos, a.End(), intStrand(a.Strand), string(stripChr(b.Chrom)), b.Pos, b.End(), intStrand(b.Strand), len(tags) - 1}
			ex.WriteAlt(&b)
		} else {
			b := bedPE{string(stripChr(b.Chrom)), b.Pos, b.End(), intStrand(b.Strand), string(stripChr(a.Chrom)), a.Pos, a.End(), intStrand(a.Strand), len(tags) - 1}
			ex.WriteAlt(&b)

		}
	}
	return len(tags)
}

func discordantByDistance(r *sam.Record, discordantDistance int) bool {
	d := iabs(r.TempLen)
	return d > discordantDistance
}

func recEnd(r *sam.Record, f *faidx.Faidx) int {
	if r.Cigar[len(r.Cigar)-1].Type() == sam.CigarSoftClipped {
		if nm, ok := r.Tag([]byte{'N', 'M'}); ok {
			v := nm[3]
			if v == 0 {
				return r.End()
			}
			/*
				seq := r.Seq.Expand()
				l := len(seq)
				chrom := r.Ref.Name()
				e := r.End()
				o := r.Cigar[len(r.Cigar)-1].Len()
				log.Println(r.Start(), r.End(), o, bigly.RefPieces(r.Pos, r.Cigar), bigly.ReadPieces(r.Cigar), r.Cigar, len(seq))
				for i := 1; i < l-1-o; i++ {
					refbase, err := f.At(chrom, e-i)
					log.Println(refbase, seq[l-i-o])
					if err != nil {
						log.Fatalf("error accessing %s:%d [%s]", chrom, e-i, err)
					}

				}
			*/
			return r.End() - int(v)
		}
	}
	return r.End()
}

// writeReferenceCoverage updates the coverage arrays (mmap'ed) for the fragment that that
// current record represents. "reference" is in the context of SVs. Reads without splitters
// or softclips are evidence for reference. The insert between reads with an insert size within
// 2SDs of the mean is also evidence for reference.
// If the pairs are discordant, it returns false if the interval is shorter than minAlignSize
func writeReferenceCoverage(r *sam.Record, fasta *faidx.Faidx, ex *excord, discordantDistance int) bool {
	ses := bigly.RefPieces(r.Pos, r.Cigar)
	n := 0
	for i := 0; i < len(ses); i += 2 {
		n += ses[i+1] - ses[i]
	}
	if n < minAlignSize {
		return false
	}

	// write coverage for the read itself here.
	// any mapped base in the read is evidence of reference at that base.
	for i := 0; i < len(ses); i += 2 {
		s, e := ses[i], ses[i+1]
		// iterate backward to remove boundschecks
		ex.updateReadCoverage(s, e)
	}
	if r.Flags&sam.Paired != sam.Paired || r.Ref.ID() != r.MateRef.ID() || r.Flags&sam.Secondary == sam.Secondary {
		return true
	}

	if discordantByDistance(r, discordantDistance) {
		return true
	}

	// pair coverage we fill from the start of left read.
	if r.Start() > r.MatePos {
		s, e := r.MatePos, recEnd(r, fasta)
		// iterate backward to remove boundschecks
		ex.updatePairCoverage(s, e)
	}
	return true
}

func stdinMain(cli *cliarg) int {
	if cli.Prefix != "" {
		panic("excord: cant specify prefix without region")
	}
	if !xopen.IsStdin() {
		panic("excord: bam streamed to sdin when no region is specified")
	}
	if cli.DiscordantDistance == 0 {
		panic("excord: when reading from stdin, you must provide a discordant distance")
	}

	ex := &excord{f: bufio.NewWriter(os.Stdout)}
	defer ex.Close()

	m := make(map[string][]*bigly.SA, 1e5)
	br, err := bam.NewReader(bufio.NewReader(os.Stdin), 2)
	pcheck(err)
	for {
		b, err := br.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			panic(err)
		}
		if uint16(b.Flags)&cli.ExcludeFlag != 0 {
			continue
		}
		if b.MapQ < cli.MinMappingQuality {
			continue
		}
		writeSplitter(b, ex)
		writeDiscordant(b, ex, cli, m)
	}

	return 0
}

func main() {
	cli := &cliarg{ExcludeFlag: uint16(sam.Unmapped | sam.QCFail | sam.Duplicate),
		MinMappingQuality: 1}
	arg.MustParse(cli)
	log.Println(cli.Region, cli.BamPath)

	if cli.Region == "" {
		os.Exit(stdinMain(cli))
	}

	chromse := strings.Split(cli.Region, ":")
	b, err := bamat.New(cli.BamPath)
	pcheck(err)
	chrom := chromse[0]
	var cLen int
	if r, ok := b.Refs[chrom]; ok {
		cLen = r.Len()
	} else if r, ok := b.Refs["chr"+chrom]; ok {
		cLen = r.Len()
	}
	if cLen <= 0 {
		panic(fmt.Sprintf("didn't find chromosome: %s\n", chrom))
	}

	if cli.DiscordantDistance == 0 {
		stats := covmed.BamInsertSizes(b.Reader, 5e5)
		b.Reader.Omit(0)
		cli.DiscordantDistance = int(stats.TemplateMean + 5*stats.TemplateSD)
		log.Printf("using distance: %d", cli.DiscordantDistance)
	}

	if cli.Prefix != "" && !strings.HasSuffix(cli.Prefix, "/") && !strings.HasSuffix(cli.Prefix, ".") {
		cli.Prefix += "."
	}

	var ex *excord

	if cli.Prefix != "" {
		ex = newExcord(cLen, cli.Prefix, cli.DiscordantDistance)
		ex.chrom = sstripChr(chrom)
	} else {
		ex = &excord{f: bufio.NewWriter(os.Stdout)}
	}
	defer ex.Close()

	var fasta *faidx.Faidx

	if cli.Fasta != "" {
		fasta, err = faidx.New(cli.Fasta)
		pcheck(err)
	}

	var start, end int
	if len(chromse) > 1 {
		se := strings.Split(chromse[1], "-")

		start, err = strconv.Atoi(se[0])
		pcheck(err)
		end, err = strconv.Atoi(se[1])
		pcheck(err)
	} else {
		start = 1
		end = cLen
	}

	it, err := b.Query(chromse[0], start-1, end)
	if err != nil { // strip chr
		if strings.HasPrefix(chromse[0], "chr") {
			it, err = b.Query(chromse[0][3:], start-1, end)
		} else {
			it, err = b.Query("chr"+chromse[0], start-1, end)
		}
	}
	pcheck(err)

	m := make(map[string][]*bigly.SA, 1e5)

	for it.Next() {
		b := it.Record()
		if uint16(b.Flags)&cli.ExcludeFlag != 0 {
			continue
		}
		if b.MapQ < cli.MinMappingQuality {
			continue
		}
		writeSplitter(b, ex)
		refs := writeDiscordant(b, ex, cli, m)
		if refs != nil && ex != nil {
			writeRefIntervals(refs, ex)
		}
		if ex != nil {
			writeReferenceCoverage(b, fasta, ex, cli.DiscordantDistance)
		}
	}
	pcheck(it.Error())
	close(ex.ch)
	ex.wg.Wait()
	s := 0
	for _, v := range ex.altMask {
		if v {
			s++
		}
	}
	log.Printf("bases covered by alts: %d, out of: %d -> %.4f%%", s, cLen, 100*float64(s)/float64(cLen))
	ex.quantize()
}

func stripChr(chrom []byte) []byte {
	if chrom[0] != 'c' {
		return chrom
	}
	if bytes.HasPrefix(chrom, []byte("chr")) {
		return chrom[3:]
	}
	return chrom
}
func sstripChr(chrom string) string {
	if chrom[0] != 'c' {
		return chrom
	}
	if strings.HasPrefix(chrom, "chr") {
		return chrom[3:]
	}
	return chrom
}
