import Foundation
import Algorithms
import Tqdm

//-------------------------------------------------------------
// Extensions, Constants, and Globals
//-------------------------------------------------------------
typealias u8 = UInt8
typealias u16 = UInt16
typealias u32 = UInt32
typealias f32 = Float 
typealias f64 = Double // lazy Rust-like abbreviations

typealias Byte = UInt8
typealias Selection = Byte
typealias Choice = Byte
typealias DieVal = Byte
typealias Slot = Byte

extension Sequence where Element: AdditiveArithmetic { func sum() -> Element { reduce(.zero, +) } }
extension String: Error {} // enables usage like: throw "error message"

let ACES:u8 = 1; let TWOS:u8 = 2; let THREES:u8 = 3;
let FOURS:u8 = 4; let FIVES:u8 = 5; let SIXES:u8 = 6;
let THREE_OF_A_KIND:u8 = 7; let FOUR_OF_A_KIND:u8 = 8;
let FULL_HOUSE:u8 = 9; let SM_STRAIGHT:u8 = 10; 
let LG_STRAIGHT:u8 = 11;
let YAHTZEE:u8 = 12; let CHANCE:u8 = 13 ;

// make an easy exponent operator
infix operator ** : ExponentiationPrecedence
precedencegroup ExponentiationPrecedence { associativity: right higherThan: MultiplicationPrecedence }
func ** (_ base: Double, _ exp: Double) -> Double { pow(base, exp) }
func ** (_ base: Float, _ exp: Float) -> Float { pow(base, exp) }
func ** (_ base: UInt, _ exp: UInt) -> UInt { UInt(pow(Float(base), Float(exp))) }
func ** (_ base: Int, _ exp: Int) -> Int { Int(pow(Float(base), Float(exp))) }


var cores = ProcessInfo().activeProcessorCount
var OUTCOME_EVS_BUFFER = (0...cores).map { _ in [f32](unsafeUninitializedCapacity:1683){$1=1683} }//  = new f32[1683,Environment.ProcessorCount];
var NEWVALS_DATA_BUFFER = (0...cores).map { _ in [u16](unsafeUninitializedCapacity:1683){$1=1683} }  //[1683,Environment.ProcessorCount]; 
var EVS_TIMES_ARRANGEMENTS_BUFFER = (0...cores).map { _ in [f32](unsafeUninitializedCapacity:1683){$1=1683} }//= new f32[1683,Environment.ProcessorCount]; 
var SORTED_DIEVALS = [DieValsID](unsafeUninitializedCapacity:32767){$1=32767}  //new DieValsID[32767];
var RANGE_IDX_FOR_SELECTION = [0,1,2,3,7,4,8,11,17,5,9,12,20,14,18,23,27,6,10,13,19,15,21,24,28,16,22,25,29,26,30,31] 
var SELECTION_RANGES = [Range<Int>](unsafeUninitializedCapacity:32){$1=32}  //new Range[32];  
var OUTCOMES = [Outcome](unsafeUninitializedCapacity:1683){$1=1683}  //new Outcome[1683]    
var OUTCOME_DIEVALS_DATA = [u16](unsafeUninitializedCapacity:1683){$1=1683} //new u16[1683]  //these 3 arrays mirror that in OUTCOMES but are contiguous and faster to access
var OUTCOME_MASK_DATA = [u16](unsafeUninitializedCapacity:1683){$1=1683}  // new u16[1683] 
var OUTCOME_ARRANGEMENTS = [f32](unsafeUninitializedCapacity:1683){$1=1683} //new f32[1683] 

var ev_cache = [ChoiceEV](unsafeUninitializedCapacity:2**30){$1=2**30}// 2^30 slots hold all unique game states 
var bar = Tqdm()

//-------------------------------------------------------------
// TEST CODE 
//-------------------------------------------------------------

// var ps = powerset([0,1,2,3,4]) 
// print(ps)

//-------------------------------------------------------------
// MAIN 
//-------------------------------------------------------------

main()

func main(){

    // TODO enable rich object formatting for debugging https://lldb.llvm.org/use/variable.html        

    //load caches
    cache_selection_ranges(); 
    cache_sorted_dievals(); 
    cache_roll_outcomes_data(); 

    // print_state_choices_header();
    let game = GameState( 
        DieVals(0,0,0,0,0), // five unrolled dice
        // DieVals(3,4,4,6,6), // five unrolled dice
        Slots(1,2,3,4,5,6,7,8,9,10,11,12,13), // all slots in an empty scorecard
        // Slots(6,12), 
        0, // current upper section total
        3, // rolls remaining
        false // yahtzee bonus available? 
    ) 

    init_bar_for(game)
    build_cache(game)

    print(ev_cache[Int(game.id)])  // # starting game state, should have expected value of 254.59
}

//-------------------------------------------------------------
//    BUILD_CACHE
//-------------------------------------------------------------

// gather up expected values in a multithreaded bottom-up fashion. this is like.. the main thing
func build_cache(_ game :GameState) {

    let range = outcomes_range_for(0b11111)
    let all_dieval_combos =  OUTCOMES[range].map {$0.dievals} 
    let placeholder_dievals = DieVals()
    let placeholder_dievals_vec = [placeholder_dievals]

    let false_true = [true, false] // NOTE These were stack alloc(?) tuples in Julia
    let just_false = [false]

    // first handle special case of the most leafy leaf calcs -- where there's one slot left and no rolls remaining
    for single_slot in game.open_slots {
        let slot = Slots([single_slot]) // set of a single slot 
        let joker_rules_in_play = (single_slot != YAHTZEE) // joker rules in effect when the yahtzee slot is not open 
        for yahtzee_bonus_available in (joker_rules_in_play ? false_true : just_false){ // yahtzee bonus -might- be available when joker rules are in play 
            for upper_total in slot.useful_upper_totals() {
                for outcome_combo in all_dieval_combos{
                    let state = GameState(outcome_combo, slot, upper_total, 0, yahtzee_bonus_available)
                    let score = state.score_first_slot_in_context()
                    let choice_ev = ChoiceEV(single_slot, f32(score))
                    ev_cache[Int(state.id)] = choice_ev
                    output(state:state, choice_ev:choice_ev)
    } } } } 

    // for each length 
    for slots_len in 1...game.open_slots.count {//Range(1, game.open_slots.Count())  {

        // for each slotset (of above length)
        for slots_vec in game.open_slots.combinations(ofCount:slots_len) {
            let slots = Slots(slots_vec);
            let joker_rules_in_play = !slots.has(YAHTZEE); // joker rules are in effect whenever the yahtzee slot is already filled 

            // for each upper total 
            for upper_total in slots.useful_upper_totals(){

                // for each yahtzee bonus possibility 
                for yahtzee_bonus_available in joker_rules_in_play ? false_true : just_false { // bonus always unavailable unless yahtzees are wild first

                    bar.update() // advance the progress bar 

                    // for each rolls remaining
                    for rolls_remaining in u8(0)...u8(3) {

                        let dieval_combos = (rolls_remaining==3 ? placeholder_dievals_vec : all_dieval_combos);

                        for dieval_combo in dieval_combos { /*Threads.@threads :static*/ 
                            process_dieval_combo(
                                rolls_remaining,
                                slots_len,
                                slots,
                                dieval_combo,
                                joker_rules_in_play,
                                yahtzee_bonus_available,
                                upper_total,
                                placeholder_dievals
                            )
    } } } } } }

}

func process_dieval_combo(_ rolls_remaining :u8, _ slots_len :Int, _ slots :Slots, _ dieval_combo :DieVals, _ joker_rules_in_play :Bool, 
                        _ yahtzee_bonus_available :Bool, _ upper_total :u8, _ placeholder_dievals :DieVals) { 

    let threadid = 0 //TODO implement actual threading // threadid = Threads.threadid()

    if (rolls_remaining==0 && slots_len > 1) { // slot selection, but not leaf calcs already recorded

        //= HANDLE SLOT SELECTION  =//

        var slot_choice_ev=ChoiceEV(0,0)

        for slot in slots { 

            // joker rules say extra yahtzees must be played in their matching upper slot if it's available
            let first_dieval = dieval_combo[1]
            let joker_rules_matter = (joker_rules_in_play && Score.yahtzee(dieval_combo)>0 && slots.has(first_dieval))
            let head_slot = (joker_rules_matter ? first_dieval : slot)

            var yahtzee_bonus_avail_now = yahtzee_bonus_available
            var upper_total_now = upper_total
            var dievals_or_placeholder = dieval_combo
            var head_plus_tail_ev:f32 = 0.0
            var rolls_remaining_now:u8 = 0

            // find the collective ev for the all the slots with this iteration's slot being first 
            // do this by summing the ev for the first (head) slot with the ev value that we look up for the remaining (tail) slots
            let passes = (slots_len==1 ? 1 : 2) // to do this, we need two passes unless there's only 1 slot left
            for i in 1...passes {
                let slots_piece = i==1 ? Slots([head_slot]) : slots.removing(head_slot)  // work on the head only or the set of slots without the head
                let upper_total_to_save = (upper_total_now + slots_piece.best_upper_total() >= 63) ? upper_total_now : (u8)(0)  // only relevant totals are cached
                let state_to_get = GameState(
                    dievals_or_placeholder,
                    slots_piece, 
                    upper_total_to_save,
                    rolls_remaining_now, 
                    yahtzee_bonus_avail_now
                )
                let choice_ev = ev_cache[Int(state_to_get.id)]
                if (i==1 && slots_len > 1) {// prep 2nd pass on relevant 1st pass only..  
                    // going into tail slots next, we may need to adjust the state based on the head choice
                    if (choice_ev.choice <= SIXES){  // adjust upper total for the next pass 
                        let added = u8(choice_ev.ev.truncatingRemainder(dividingBy: 100.0) ) // the modulo 100 here removes any yathzee bonus from ev since that doesnt' count toward upper bonus total
                        upper_total_now = u8(min(63, upper_total_now + added))
                    } else if (choice_ev.choice==YAHTZEE) {  // adjust yahtzee related state for the next pass
                        if (choice_ev.ev>0.0) {yahtzee_bonus_avail_now=true}
                    } 
                    rolls_remaining_now=3 // for upcoming tail lookup, we always want the ev for 3 rolls remaining
                    dievals_or_placeholder = placeholder_dievals // for 4 rolls remaining, use "wildcard" representative dievals since dice don't matter when rolling all of them
                }
                head_plus_tail_ev += choice_ev.ev
            }//for i in passes 

            if head_plus_tail_ev >= slot_choice_ev.ev { 
                slot_choice_ev = ChoiceEV(slot, head_plus_tail_ev)
            } 
            
            if joker_rules_matter {break} // if joker-rules-matter we were forced to choose one slot, so we can skip trying the rest  

        }//for slot in slots                               

        let state_to_set = GameState(
            dieval_combo,
            slots,
            upper_total, 
            0, 
            yahtzee_bonus_available
        ); 
        ev_cache[Int(state_to_set.id)] = slot_choice_ev;
        // output_state_choice(this, state_to_set, slot_choice_ev)

    } else if (rolls_remaining > 0) {  

    //= HANDLE DICE SELECTION =//    

        let next_roll = rolls_remaining-1
        var best = ChoiceEV(0,0.0)//  selections are bitfields where '1' means roll and '0' means don't roll 
        let selections = rolls_remaining==3 ? u8(0b11111)...u8(0b11111) : u8(0b00000)...u8(0b11111)//TODO test hoist? //select all dice on the initial roll, else try all selections
        
        for selection in selections { // we'll try each selection against this starting dice combo  
            let avg_ev_for_selection = avg_ev(dieval_combo, selection, slots, upper_total, next_roll,yahtzee_bonus_available, ev_cache, threadid) // @inline
            if (avg_ev_for_selection > best.ev){
                best = ChoiceEV(selection, avg_ev_for_selection)
            } 
        } 
        let state_to_set = GameState(
            dieval_combo,
            slots, 
            upper_total, 
            rolls_remaining, 
            yahtzee_bonus_available 
        ) 
        output(state:state_to_set, choice_ev:best);
        ev_cache[Int(state_to_set.id)] = best;

    }// if rolls_remaining...  
}// process_dieval_combo


func avg_ev(_ start_dievals :DieVals, _ selection :Selection, _ slots :Slots, _ upper_total :u8, 
            _ next_roll :u8, _ yahtzee_bonus_available :Bool, _ ev_cache :[ChoiceEV] , _ threadid:Int) -> f32 { 

    var total_ev_for_selection:f32 = 0.0 ;
    var outcomes_arrangements_count:f32 = 0.0;
    let range = outcomes_range_for(selection);

    for i in range { //@inbounds @simd
        NEWVALS_DATA_BUFFER[threadid][i] = (u16)(start_dievals.data & OUTCOME_MASK_DATA[i]);
        NEWVALS_DATA_BUFFER[threadid][i] |= OUTCOME_DIEVALS_DATA[i];
    } 

    let floor_state_id = GameState(
        DieVals(),
        slots, 
        upper_total, 
        next_roll, // we'll average all the 'next roll' possibilities (which we'd calclated last) to get ev for 'this roll' 
        yahtzee_bonus_available 
    ).id

    for i in range { //@inbounds @simd
        //= gather sorted =#
            let u = i
            let newvals_datum = Int(NEWVALS_DATA_BUFFER[threadid][u])
            let sorted_dievals_id  = SORTED_DIEVALS[newvals_datum].id 
        //= gather ev =#
            let state_to_get_id = Int(floor_state_id) | Int(sorted_dievals_id)
            let cache_entry = ev_cache[state_to_get_id]
            OUTCOME_EVS_BUFFER[threadid][u] = cache_entry.ev
    } 

    for i in range { //@inbounds @simd
    // foreach(int i in range) {// we looped through die "combos" but we need to average all "perumtations" // @fastmath @inbounds @simd ivdep 
        EVS_TIMES_ARRANGEMENTS_BUFFER[threadid][i] = OUTCOME_EVS_BUFFER[threadid][i] * OUTCOME_ARRANGEMENTS[i]
        total_ev_for_selection +=  EVS_TIMES_ARRANGEMENTS_BUFFER[threadid][i]
        outcomes_arrangements_count += OUTCOME_ARRANGEMENTS[i] 
    } 

    return  total_ev_for_selection / outcomes_arrangements_count;

} // avg_ev

//-------------------------------------------------------------
// INITIALIZERS 
//-------------------------------------------------------------

// this generates the ranges that correspond to the outcomes, within the set of all outcomes, indexed by a give selection 
func cache_selection_ranges() {

    var s = 0
    SELECTION_RANGES[0] = 0..<1
    let combos = powerset([0,1,2,3,4])// (new List<int>(){0,1,2,3,4}).powerset();

    var i = 0
    for combo in combos {
        let count = Int( n_take_r(6, UInt(combo.count), order_matters:false, with_replacement:true) )
        SELECTION_RANGES[i] = s..<(s + count)
        s += count
        i+=1
    } 
}

// for fast access later, this generates an array of dievals in sorted form, 
// along with each's unique "ID" between 0-252, indexed by DieVals.data
func cache_sorted_dievals() { 
    
    // TODO move all this inside block above?
    SORTED_DIEVALS[0] = DieValsID(); // first one is for the special wildcard 
    let one_to_six = [u8](1...6)
    let combos = combos_with_rep(one_to_six,5) 
    for (i, combo) in combos.enumerated() {
        let dv_combo = DieVals(combo);
        let int_perms = combo.uniquePermutations()
        let dv_perms = int_perms.map { DieVals($0) } 
        for perm in dv_perms { SORTED_DIEVALS[Int(perm.data)] = DieValsID(dv_combo, u8(i)) }
    }
}

//preps the caches of roll outcomes data for every possible 5-die selection, where '0' represents an unselected die """
func cache_roll_outcomes_data() { 

    var i=0 
    let idx_combos = powerset([Int](0...4))
    let one_thru_six:[u8] = [1,2,3,4,5,6]
    for idx_combo_vec in idx_combos { 
        var dievals_vec:[DieVal] = [0,0,0,0,0] // new DieVal[5]
        let die_count = idx_combo_vec.count
        for dievals_combo_vec in combos_with_rep(one_thru_six, die_count) {
            var mask_vec:[u8] = [0b111,0b111,0b111,0b111,0b111]
            for (j, val) in dievals_combo_vec.enumerated() {
                let idx = idx_combo_vec[j]
                dievals_vec[idx] = DieVal(val) 
                mask_vec[idx]=DieVal()
            } 
            let arrangement_count = distinct_arrangements_for(dievals_combo_vec)
            let dievals = DieVals(dievals_vec)
            let mask = DieVals(mask_vec)
            OUTCOME_DIEVALS_DATA[i] = dievals.data
            OUTCOME_MASK_DATA[i] = mask.data
            OUTCOME_ARRANGEMENTS[i] = arrangement_count
            OUTCOMES[i] = Outcome( dievals, mask, arrangement_count)
            i+=1;
        } 
    } 
} 

func init_bar_for(_ game :GameState) {
    bar = Tqdm(total:game.counts())
} 

func output(state :GameState, choice_ev :ChoiceEV ){ 
    // Uncomment below for more verbose progress output at the expense of speed 
    // print_state_choice(state, choice_ev)
} 


//-------------------------------------------------------------
//  UTILS
//-------------------------------------------------------------

func powerset<T>(_ set :[T] ) -> [[T]]{
    let size = UInt(set.count)
    let setsize :UInt = 2**size-1//set_size of power set of a set with set_size n is (2**n -1)
    var outerList = [[T]]()  
    for i in 0...setsize {// Run from counter 000..0 to 111..1
        var innerList = [T]() // Check each jth bit in the counter is set If set then add jth element from set 
        for j in 0..<size {
            if ((i & (1 << j)) > 0) {innerList.append(set[Int(j)])}
        }
        outerList.append(innerList)
    }
    return outerList;
}

// count of arrangements that can be formed from r selections, chosen from n items, 
// where order DOES or DOESNT matter, and WITH or WITHOUT replacement, as specified.
func n_take_r(_ n:UInt , _ r:UInt, order_matters:Bool=false, with_replacement:Bool=false) -> UInt{  
    if (order_matters){  //  order matters; we're counting "permutations" 
        if (with_replacement) {
            return n**r;
        } else { //  no replacement
            return factorial(n) / factorial(n-r);  //  this = factorial(n) when r=n
        }
    } else { //  we're counting "combinations" where order doesn't matter; there are less of these 
        if (with_replacement) {
            return factorial(n+r-1) / (factorial(r)*factorial(n-1));
        } else { //  no replacement
            return factorial(n) / (factorial(r)*factorial(n-r));
        }
    }
}

func factorial(_ n:UInt) -> UInt { 
    if (n<=1) {return 1}
    return n * factorial(n-1);
}

func combos_with_rep<T>(elements: ArraySlice<T>, k: Int) -> [[T]] {
    if k == 0 { return [[]] }
    guard let first = elements.first else { return [] }
    let head = [first]
    let subcombos = combos_with_rep(elements: elements, k: k - 1)
    var ret = subcombos.map { head + $0 }
    ret += combos_with_rep(elements: elements.dropFirst(), k: k)
    return ret
}

func combos_with_rep<T>(_ elements: Array<T>, _ k: Int) -> [[T]] {
    return combos_with_rep(elements: ArraySlice(elements), k: k)
}

func distinct_arrangements_for(_ dieval_vec:[DieVal]) -> f32 { 
    // var key_counts = dieval_vec.GroupBy(x=>x).Select(g=>(g.Key, (u8)g.Count()));
    let key_counts = Dictionary(grouping: dieval_vec, by: { $0 }).mapValues { $0.count } // group val with counts
    var divisor:UInt=1
    var non_zero_dievals:UInt=0
    for (key, count) in key_counts{  
        if (key != 0){  
            divisor *= factorial(UInt(count))
            non_zero_dievals += UInt(count)
        } 
    } 
    return f32( factorial(non_zero_dievals) / divisor )
} 

// returns a range which corresponds the precomputed dice roll outcome data corresponding to the given selection
func outcomes_range_for(_ selection :Selection) -> Range<Int>{
    let idx = RANGE_IDX_FOR_SELECTION[Int(selection)];
    let range = SELECTION_RANGES[idx]; // for @inbounds, see https://blog.tedd.no/2020/06/01/faster-c-array-access/
    return range;
} 

func print_state_choices_header() { 
    print("choice_type,choice,dice,rolls_remaining,upper_total,yahtzee_bonus_avail,open_slots,expected_value");
} 

// should produce one line of output kinda like ...
// D,01111,65411,2,31,Y,1_3_4_6_7_8_11_,119.23471
// S,13,66641,0,11,Y,3_4_5_9_10_13_,113.45208
func print_state_choice(_ state :GameState , _ choice_ev: ChoiceEV ) { 
    let Y="Y"; let N="N"; let S="S"; let D="D"; let C=","; // TODO hoist these to constants
    var sb:String=""; sb.reserveCapacity(60)
    if (state.rolls_remaining==0){ // for slot choice 
        sb += (S); sb+=(C);
        sb += (choice_ev.choice.description); // for dice choice 
    } else {
        sb+=(D); sb+=(C);
        sb+=("00000"+choice_ev.choice.description).suffix(5)
    }
    sb+=(C);
    sb+=(state.sorted_dievals.description); sb+=(C);
    sb+=(state.rolls_remaining.description); sb+=(C);
    sb+=(state.upper_total.description); sb+=(C);
    sb+=(state.yahtzee_bonus_avail ? Y : N); sb+=(C);
    sb+=(state.open_slots.description); sb+=(C);
    sb+=(choice_ev.ev.description);
    // sb+=(C); sb+=(state.id.description);
    print(sb);
} 


//-------------------------------------------------------------
//DieValsID
//-------------------------------------------------------------
struct DieValsID { 
    var dievals = DieVals() 
    var id:u8
    init(_ dievals:DieVals, _ id:u8 ){ self.dievals = dievals; self.id = id; }
    init() { dievals = DieVals(); self.id = 0; }
}

//-------------------------------------------------------------
// Outcome
//-------------------------------------------------------------
struct Outcome { 
    var dievals:DieVals = DieVals()
    var mask:DieVals = DieVals()// stores a pre-made mask for blitting this outcome onto a GameState.DieVals.data u16 later
    var arrangements:f32 = 0.0 // how many indistinguisable ways can these dievals be arranged (ie swapping identical dievals)
    init() { } 
    init(_ dievals:DieVals, _ mask:DieVals, _ arrangements:f32) { 
        self.dievals = dievals 
        self.mask = mask 
        self.arrangements = arrangements 
    }
}

//#=-------------------------------------------------------------
//DieVals
//-------------------------------------------------------------=#

// the following LLDB command will format DieVals with meaningful values in the debugger 
/*    
type summary add swiftzbot.DieVals --python-script 
"v = valobj.GetChildMemberWithName('data').GetValueAsUnsigned(0); 
d1=str(v & 0b_111); d2=str((v & 0b_111_000)>>3); d3=str((v & 0b_111_000_000)>>6); 
d4=str((v & 0b_111_000_000_000)>>9); d5=str((v & 0b_111_000_000_000_000)>>12); return d1+d2+d3+d4+d5" 
*/

struct DieVals : Collection, CustomStringConvertible, CustomDebugStringConvertible  
    { // TODO more performant to make this lightweight Julia-like struct without the iterator baggage?
    public var data :u16 = 0// 5 dievals, each from 0 to 6, can be encoded in 2 bytes total, each taking 3 bits

    //initializers
    init() {}
    init(_ args:DieVal... ) { self.init(args) }
    init(_ dievals :[DieVal] ) {
        for (i, dieval) in dievals.enumerated() { 
            self.data |= u16(dieval) << (i*3)  
        } 
    }

    //implement Collection protocol
    typealias Index = Int
    typealias Element = DieVal

    var startIndex: Int { return 0 }
    var endIndex: Int { return 5 }

    subscript(position: Int) -> DieVal {
        get { return DieVal( (self.data >> (position*3)) & 0b111) }
    }

    func index(after i: Int) -> Int { return i + 1 }

    var debugDescription: String { return description }

    //implement CustomStringConvertible protocol
    public var description: String { return "\(self[4])\(self[3])\(self[2])\(self[1])\(self[0])" }

}

//-------------------------------------------------------------
// SLOTS 
// ------------------------------------------------------------
// the following LLDB command will format SortedSlots with meaningful values in the debugger 
// type summary add swiftzbot.Slots --summary-string "${var.data%b}" 

struct Slots : Collection, CustomStringConvertible, CustomDebugStringConvertible  
    {

    public var data:u16 = 0 // 13 sorted Slots can be positionally encoded in one u16

    init() {}

    init(_ args :Slot... ) {  self.init(args) }

    init(_ array :[Slot]){  
        var mask :u16 
        for slot in array {
            mask = u16(0x0001) << slot
            self.data |= mask; // force on
        } 
    }

    //implement Collection protocol
    typealias Index = Int
    typealias Element = Slot

    var startIndex: Int { return 0 }
    var endIndex: Int { return data.nonzeroBitCount}

    subscript(position: Int) -> Slot {
        get { 
            // @assert(i<=length(self))
            var bits = data
            var bit_index=0
            let i = position+1
            var j=1 //the slots.data does not use the rightmost (0th) bit 
            while (j <= i){ 
                bit_index = bits.trailingZeroBitCount
                bits = (bits & (~( 1 << bit_index) ))  //unset bit
                j+=1
            } 
            return Slot(bit_index)
        }
    }

    func index(after i: Int) -> Int { 
        return i + 1 
    }

    public func has(_ slot:Slot) -> Bool { 
        return 0x0000 < (data & (0x0001<<(u16)(slot)));
    } 

    public func removing(_ slot_to_remove :Slot) -> Slots  { 
        let mask:u16 = ~( 1 << u16(slot_to_remove) );
        let newdata:u16 = u16(self.data & mask); //# force off
        return slotsFromData(newdata);
    } 

    private func slotsFromData(_ data:u16) -> Slots {
        var slots = Slots(); 
        slots.data = data;
        return slots;
    }

    private func iseven(_ x :Slot)->Bool { x%2==0 } 
    private func iseven(_ x :Int)->Bool { x%2==0 } // TODO some generic way to avoid dupe ?

    // a non-exact but fast estimate of relevant_upper_totals
    // ie the unique and relevant "upper bonus total" that could have occurred from the previously used upper slots
    func useful_upper_totals() -> [u8] { 
        var totals:[u8] = Array(0...63) //(x for x in 0:63)
        let used_uppers = self.used_upper_slots()
        if (used_uppers.allSatisfy(iseven))  {
            totals = totals.filter { iseven($0) } 
        } 
        // filter out the lowish totals that aren't relevant because they can't reach the goal with the upper slots remaining 
        // this filters out a lot of dead state space but means the lookup function must later map extraneous deficits to a default 
        let best_unused_slot_total = self.best_upper_total()
        // totals = (x for x in totals if x + best_unused_slot_total >=63 || x==0)
        // totals = from x in totals where (x + best_unused_slot_total >=63 || x==0) select x
        totals = totals.filter { $0 + best_unused_slot_total >= 63 || $0==0 }
        return totals;
    }

    func best_upper_total() -> u8 { 
        var sum:u8 = 0;
        for x in self { if(6<x){break}; sum+=x } 
        return sum*5
    } 

    func used_upper_slots() -> Slots {
        let all_bits_except_unused_uppers = ~self.data; // "unused" slots (as encoded in .data) are not "previously used", so blank those out
        let all_upper_slot_bits = (u16) ((1<<7)-2);  // upper slot bits are those from 2^1 through 2^6 (.data encoding doesn't use 2^0)
        let previously_used_upper_slot_bits = (u16) (all_bits_except_unused_uppers & all_upper_slot_bits);
        return slotsFromData( previously_used_upper_slot_bits );
    } 

    func best_upper_total(slots:Slots) -> u8 {  
        var sum:u8=0
        for s in slots { if (6<s) {break}; sum+=s } 
        return u8(sum*5)
    } 

    //implement CustomStringConvertible protocol
    public var description: String { 
        var sb:String=""; sb.reserveCapacity(13)
        for s in self { sb += "\(s)_" }
        return sb
    }

    var debugDescription: String { return description }

}

//-------------------------------------------------------------
//ChoiceEV
//-------------------------------------------------------------
struct ChoiceEV : CustomStringConvertible {
    public var choice :Choice
    public var ev :f32 ;
    public var description: String { return "ChoiceEV(\(choice),\(ev))" }
    init(_ choice:Choice, _ ev:f32) { self.choice=choice; self.ev=ev }
}

//-------------------------------------------------------------
//GameState
//-------------------------------------------------------------

struct GameState {
    public var id:u32 = 0; // 30 bits # with the id, 
    //we can store all of below in a sparse array using 2^(8+13+6+2+1) entries = 1_073_741_824 entries = 5.2GB when storing 40bit ResultEVs 
    public var sorted_dievals :DieVals ;// 15 bits OR 8 bits once convereted to a DieValID (252 possibilities)
    public var open_slots :Slots ;// 13 bits        = 8_192 possibilities
    public var upper_total :u8 ;// = 6 bits         = 64    "
    public var rolls_remaining :u8 ;// = 2 bits     = 4     "  
    public var yahtzee_bonus_avail :Bool ;// = 1bit = 2     "

    init(_ sorted_dievals:DieVals, _ open_slots:Slots, _ upper_total:u8, _ rolls_remaining:u8 , _ yahtzee_bonus_avail:Bool) { 
        let dievals_id:u8  = SORTED_DIEVALS[Int(sorted_dievals.data)].id // this is the 8-bit encoding of self.sorted_dievals
        self.id = u32(dievals_id)                  // self.id will use 30 bits total...
        self.id |= ( u32(open_slots.data)        << 7)   // slots.data doesn't use its rightmost bit so we only shift 7 to make room for the 8-bit dieval_id above 
        self.id |= ( u32(upper_total)            << 21)  // make room for 13+8 bit stuff above 
        self.id |= ( u32(rolls_remaining)        << 27)  // make room for the 13+8+6 bit stuff above
        self.id |= ( u32(yahtzee_bonus_avail ? 1 : 0) << 29)   // make room for the 13+8+6+2 bit stuff above
        self.sorted_dievals = sorted_dievals
        self.open_slots = open_slots
        self.upper_total = upper_total
        self.rolls_remaining = rolls_remaining
        self.yahtzee_bonus_avail = yahtzee_bonus_avail
    } 

    // calculate relevant counts for gamestate: required lookups and saves
    func counts() -> Int { 
        var ticks = 0 
        let false_true = [true, false]
        let just_false = [false]
        for subset_len in 1...open_slots.count { //Range(1,open_slots.Count)){
            let combos = open_slots.combinations(ofCount: subset_len)
            for slots_vec in combos { 
                let slots = Slots(slots_vec)
                let joker_rules = slots.has(YAHTZEE) // yahtzees aren't wild whenever yahtzee slot is still available 
                let totals = slots.useful_upper_totals() 
                for _ in totals {
                    for _ in joker_rules ? false_true : just_false {
                        // var slot_lookups = (subset_len * subset_len==1? 1 : 2) * 252 // * subset_len as u64
                        // var dice_lookups = 848484 // // previoiusly verified by counting up by 1s in the actual loop. however chunking forward is faster 
                        // lookups += (dice_lookups + slot_lookups) this tends to overflow so use "normalized" ticks below
                        ticks+=1 // this just counts the cost of one pass through the bar.tick call in the dice-choose section of build_cache() loop
        } } } } 
        return Int(ticks)
    } 

    public func score_first_slot_in_context() -> u8 { 

        // score slot itself w/o regard to game state 
            var it = open_slots.makeIterator()
            let slot = it.next()! // first slot in open_slots
            var score = Score.slot_with_dice(slot, sorted_dievals) 

        // add upper bonus when needed total is reached 
            if (slot<=SIXES && upper_total<63){
                let new_total = min(upper_total+score, 63 ) 
                if (new_total==63) { // we just reach bonus threshold
                    score += 35   // add the 35 bonus points 
                }
            } 

        // special handling of "joker rules" 
            let just_rolled_yahtzee = Score.yahtzee(sorted_dievals)==50
            let joker_rules_in_play = (slot != YAHTZEE) // joker rules in effect when the yahtzee slot is not open 
            if (just_rolled_yahtzee && joker_rules_in_play){ // standard scoring applies against the yahtzee dice except ... 
                if (slot==FULL_HOUSE) {score=25}
                if (slot==SM_STRAIGHT){score=30}
                if (slot==LG_STRAIGHT){score=40}
            } 

        // # special handling of "extra yahtzee" bonus per rules
            if (just_rolled_yahtzee && yahtzee_bonus_avail) {score+=100}

        return score
    } 

} 

//-------------------------------------------------------------
//SCORING FNs
//-------------------------------------------------------------
struct Score {

    public static func upperbox(_ boxnum :u8, _ sorted_dievals :DieVals) -> u8 { 
        var sum:u8=0
        for d in sorted_dievals {if (d==boxnum) {sum+=boxnum} }
        return sum 
    } 

    public static func n_of_a_kind(_ n:u8, _ sorted_dievals:DieVals) -> u8 { 
        var inarow:u8=1; var maxinarow:u8=1; var lastval:u8=100; var sum:u8=0
        for x in sorted_dievals { 
            if (x==lastval && x != 0) {inarow += 1} else {inarow=1}; 
            maxinarow = max(inarow,maxinarow);
            lastval = x;
            sum+=x;
        } 
        if (maxinarow>=n) {return sum} else {return 0} ;
    } 

    public static func straight_len(_ sorted_dievals :DieVals) -> u8 { 
        var inarow :u8=1
        var lastval :u8=254 // stub
        var maxinarow :u8=1
        for x in sorted_dievals {
            if (x==lastval+1 && x != 0){ 
                inarow+=1
            } else if (x != lastval) { 
                inarow=1 
            }
            maxinarow = max(inarow, maxinarow)
            lastval = x
        } 
        return maxinarow
    } 

    public static func aces(_ sorted_dievals :DieVals) -> u8             { return upperbox(0x1,sorted_dievals)}
    public static func twos(_ sorted_dievals :DieVals) -> u8             { return upperbox(0x2,sorted_dievals)} 
    public static func threes(_ sorted_dievals :DieVals) -> u8           { return upperbox(0x3,sorted_dievals)} 
    public static func fours(_ sorted_dievals :DieVals) -> u8            { return upperbox(0x4,sorted_dievals)} 
    public static func fives(_ sorted_dievals :DieVals) -> u8            { return upperbox(0x5,sorted_dievals)} 
    public static func sixes(_ sorted_dievals :DieVals) -> u8            { return upperbox(0x6,sorted_dievals)} 
        
    public static func three_of_a_kind(_ sorted_dievals :DieVals) -> u8  { return n_of_a_kind(0x3,sorted_dievals)} 
    public static func four_of_a_kind(_ sorted_dievals :DieVals) -> u8   { return n_of_a_kind(0x4,sorted_dievals)} 
    public static func sm_str8(_ sorted_dievals :DieVals) -> u8          { return (u8)(straight_len(sorted_dievals)>=4 ? 30 : 0)}
    public static func lg_str8(_ sorted_dievals :DieVals) -> u8          { return (u8)(straight_len(sorted_dievals)==5 ? 40 : 0)}

    // The official rule is that a Full House is "three of one number and two of another
    public static func fullhouse(_ sorted_dievals :DieVals) -> u8 { 
        var valcounts1 = 0; var valcounts2 = 0;
        var j=0;
        for (i,val) in sorted_dievals.enumerated() { 
            if (val==0) {return u8(0) }
            if (j==0 || sorted_dievals[i] != sorted_dievals[i-1]) {j+=1} 
            if (j==1) {valcounts1+=1} 
            if (j==2) {valcounts2+=1} 
            if (j==3) {return 0 }
        } 
        if (valcounts1==3 && valcounts2==2 || valcounts2==3 && valcounts1==2) {return 25} 
        return 0 
    } 
        
    public static func chance(_ sorted_dievals :DieVals) -> u8 {return u8( sorted_dievals.sum() )}
        
    public static func yahtzee(_ sorted_dievals :DieVals) -> u8 { 
        if (sorted_dievals[0]==0) {return 0} ; 
        return (u8)(sorted_dievals[0] == sorted_dievals[4] ? 50 : 0) ;
    }

    // reports the score for a set of dice in a given slot w/o regard for exogenous gamestate (bonuses, yahtzee wildcards etc) 
    public static func slot_with_dice(_ slot :Slot, _ sorted_dievals :DieVals) -> u8 { 
        if (slot==ACES) {return aces(sorted_dievals)}  
        if (slot==TWOS) {return twos(sorted_dievals)}  
        if (slot==THREES) {return threes(sorted_dievals)}  
        if (slot==FOURS) {return fours(sorted_dievals)}  
        if (slot==FIVES) {return fives(sorted_dievals)}  
        if (slot==SIXES) {return sixes(sorted_dievals)}  
        if (slot==THREE_OF_A_KIND) {return three_of_a_kind(sorted_dievals)}  
        if (slot==FOUR_OF_A_KIND) {return four_of_a_kind(sorted_dievals)}  
        if (slot==SM_STRAIGHT) {return Score.sm_str8(sorted_dievals)}  
        if (slot==LG_STRAIGHT) {return Score.lg_str8(sorted_dievals)}  
        if (slot==FULL_HOUSE) {return fullhouse(sorted_dievals)}  
        if (slot==CHANCE) {return chance(sorted_dievals)}  
        if (slot==YAHTZEE) {return yahtzee(sorted_dievals)}  
        assert(false) // shouldn't get here
        return 0
    }

}

