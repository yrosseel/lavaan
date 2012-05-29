


<!DOCTYPE html>
<html>
  <head prefix="og: http://ogp.me/ns# fb: http://ogp.me/ns/fb# githubog: http://ogp.me/ns/fb/githubog#">
    <meta charset='utf-8'>
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
        <title>lavaan/R/02lavaanUser.R at master · yrosseel/lavaan · GitHub</title>
    <link rel="search" type="application/opensearchdescription+xml" href="/opensearch.xml" title="GitHub" />
    <link rel="fluid-icon" href="https://github.com/fluidicon.png" title="GitHub" />
    <link rel="shortcut icon" href="/favicon.ico" type="image/x-icon" />

    
    

    <meta content="authenticity_token" name="csrf-param" />
<meta content="sNeiC2x0W47+BZLr5I85XeeqsfGGmIGAXtaGQKx+Q4k=" name="csrf-token" />

    <link href="https://a248.e.akamai.net/assets.github.com/stylesheets/bundles/github-e2a66f70d831aada662d99108fbb8d26c1efb6eb.css" media="screen" rel="stylesheet" type="text/css" />
    <link href="https://a248.e.akamai.net/assets.github.com/stylesheets/bundles/github2-e3d528062f9480d9819936572554df3183b4e6b2.css" media="screen" rel="stylesheet" type="text/css" />
    
    


    <script src="https://a248.e.akamai.net/assets.github.com/javascripts/bundles/frameworks-a450c7f907bdc1ee6b362ab1ecca811c761fd259.js" type="text/javascript"></script>
    
    <script defer="defer" src="https://a248.e.akamai.net/assets.github.com/javascripts/bundles/github-8d478b2120b2fe14f438e1591dd2e2ac86f44d48.js" type="text/javascript"></script>
    
    

      <link rel='permalink' href='/yrosseel/lavaan/blob/79b51b6496090dca20de417c2b60c57d3f0c93d2/R/02lavaanUser.R'>
    <meta property="og:title" content="lavaan"/>
    <meta property="og:type" content="githubog:gitrepository"/>
    <meta property="og:url" content="https://github.com/yrosseel/lavaan"/>
    <meta property="og:image" content="https://a248.e.akamai.net/assets.github.com/images/gravatars/gravatar-140.png?1334862345"/>
    <meta property="og:site_name" content="GitHub"/>
    <meta property="og:description" content="lavaan - an R package for structural equation modeling and more"/>

    <meta name="description" content="lavaan - an R package for structural equation modeling and more" />

  <link href="https://github.com/yrosseel/lavaan/commits/master.atom" rel="alternate" title="Recent Commits to lavaan:master" type="application/atom+xml" />

  </head>


  <body class="logged_out page-blob linux vis-public env-production " data-blob-contribs-enabled="yes">
    <div id="wrapper">

    
    
    

      <div id="header" class="true clearfix">
        <div class="container clearfix">
          <a class="site-logo" href="https://github.com/">
            <!--[if IE]>
            <img alt="GitHub" class="github-logo" src="https://a248.e.akamai.net/assets.github.com/images/modules/header/logov7.png?1334862346" />
            <img alt="GitHub" class="github-logo-hover" src="https://a248.e.akamai.net/assets.github.com/images/modules/header/logov7-hover.png?1334862346" />
            <![endif]-->
            <img alt="GitHub" class="github-logo-4x" height="30" src="https://a248.e.akamai.net/assets.github.com/images/modules/header/logov7@4x.png?1337118066" />
            <img alt="GitHub" class="github-logo-4x-hover" height="30" src="https://a248.e.akamai.net/assets.github.com/images/modules/header/logov7@4x-hover.png?1337118066" />
          </a>

                  <!--
      make sure to use fully qualified URLs here since this nav
      is used on error pages on other domains
    -->
    <ul class="top-nav logged_out">
        <li class="pricing"><a href="https://github.com/plans">Signup and Pricing</a></li>
        <li class="explore"><a href="https://github.com/explore">Explore GitHub</a></li>
      <li class="features"><a href="https://github.com/features">Features</a></li>
        <li class="blog"><a href="https://github.com/blog">Blog</a></li>
      <li class="login"><a href="https://github.com/login?return_to=%2Fyrosseel%2Flavaan%2Fblob%2Fmaster%2FR%2F02lavaanUser.R">Login</a></li>
    </ul>



          
        </div>
      </div>

      

            <div class="site hfeed" itemscope itemtype="http://schema.org/WebPage">
      <div class="container hentry">
        <div class="pagehead repohead instapaper_ignore readability-menu">
        <div class="title-actions-bar">
          



              <ul class="pagehead-actions">



          <li>
            <span class="watch-button"><a href="/login?return_to=%2Fyrosseel%2Flavaan" class="minibutton btn-watch js-toggler-target entice tooltipped leftwards" title="You must be logged in to use this feature" rel="nofollow"><span class="icon"></span> Watch</a><a class="social-count js-social-count" href="/yrosseel/lavaan/watchers">14</a></span>
          </li>
          <li>
            <a href="/login?return_to=%2Fyrosseel%2Flavaan" class="minibutton btn-fork js-toggler-target fork-button entice tooltipped leftwards"  title="You must be logged in to use this feature" rel="nofollow"><span class="icon"></span>Fork</a><a href="/yrosseel/lavaan/network" class="social-count">5</a>
          </li>

    </ul>

          <h1 itemscope itemtype="http://data-vocabulary.org/Breadcrumb" class="entry-title public">
            <span class="repo-label"><span>public</span></span>
            <span class="mega-icon mega-icon-public-repo"></span>
            <span class="author vcard">
<a href="/yrosseel" class="url fn" itemprop="url" rel="author">              <span itemprop="title">yrosseel</span>
              </a></span> /
            <strong><a href="/yrosseel/lavaan" class="js-current-repository">lavaan</a></strong>
          </h1>
        </div>

          

  <ul class="tabs">
    <li><a href="/yrosseel/lavaan" class="selected" highlight="repo_sourcerepo_downloadsrepo_commitsrepo_tagsrepo_branches">Code</a></li>
    <li><a href="/yrosseel/lavaan/network" highlight="repo_network">Network</a>
    <li><a href="/yrosseel/lavaan/pulls" highlight="repo_pulls">Pull Requests <span class='counter'>1</span></a></li>

      <li><a href="/yrosseel/lavaan/issues" highlight="repo_issues">Issues <span class='counter'>2</span></a></li>


    <li><a href="/yrosseel/lavaan/graphs" highlight="repo_graphsrepo_contributors">Graphs</a></li>

  </ul>
 
<div class="frame frame-center tree-finder" style="display:none"
      data-tree-list-url="/yrosseel/lavaan/tree-list/79b51b6496090dca20de417c2b60c57d3f0c93d2"
      data-blob-url-prefix="/yrosseel/lavaan/blob/79b51b6496090dca20de417c2b60c57d3f0c93d2"
    >

  <div class="breadcrumb">
    <span class="bold"><a href="/yrosseel/lavaan">lavaan</a></span> /
    <input class="tree-finder-input js-navigation-enable" type="text" name="query" autocomplete="off" spellcheck="false">
  </div>

    <div class="octotip">
      <p>
        <a href="/yrosseel/lavaan/dismiss-tree-finder-help" class="dismiss js-dismiss-tree-list-help" title="Hide this notice forever" rel="nofollow">Dismiss</a>
        <span class="bold">Octotip:</span> You've activated the <em>file finder</em>
        by pressing <span class="kbd">t</span> Start typing to filter the
        file list. Use <span class="kbd badmono">↑</span> and
        <span class="kbd badmono">↓</span> to navigate,
        <span class="kbd">enter</span> to view files.
      </p>
    </div>

  <table class="tree-browser" cellpadding="0" cellspacing="0">
    <tr class="js-header"><th>&nbsp;</th><th>name</th></tr>
    <tr class="js-no-results no-results" style="display: none">
      <th colspan="2">No matching files</th>
    </tr>
    <tbody class="js-results-list js-navigation-container">
    </tbody>
  </table>
</div>

<div id="jump-to-line" style="display:none">
  <h2>Jump to Line</h2>
  <form accept-charset="UTF-8">
    <input class="textfield" type="text">
    <div class="full-button">
      <button type="submit" class="classy">
        <span>Go</span>
      </button>
    </div>
  </form>
</div>


<div class="subnav-bar">

  <ul class="actions subnav">
    <li><a href="/yrosseel/lavaan/tags" class="blank" highlight="repo_tags">Tags <span class="counter">0</span></a></li>
    <li><a href="/yrosseel/lavaan/downloads" class="blank downloads-blank" highlight="repo_downloads">Downloads <span class="counter">0</span></a></li>
    
  </ul>

  <ul class="scope">
    <li class="switcher">

      <div class="context-menu-container js-menu-container js-context-menu">
        <a href="#"
           class="minibutton bigger switcher js-menu-target js-commitish-button btn-branch repo-tree"
           data-hotkey="w"
           data-master-branch="master"
           data-ref="master">
           <span><span class="icon"></span><i>branch:</i> master</span>
        </a>

        <div class="context-pane commitish-context js-menu-content">
          <a href="javascript:;" class="close js-menu-close"><span class="mini-icon mini-icon-remove-close"></span></a>
          <div class="context-title">Switch Branches/Tags</div>
          <div class="context-body pane-selector commitish-selector js-navigation-container">
            <div class="filterbar">
              <input type="text" id="context-commitish-filter-field" class="js-navigation-enable" placeholder="Filter branches/tags" data-filterable />

              <ul class="tabs">
                <li><a href="#" data-filter="branches" class="selected">Branches</a></li>
                <li><a href="#" data-filter="tags">Tags</a></li>
              </ul>
            </div>

            <div class="js-filter-tab js-filter-branches" data-filterable-for="context-commitish-filter-field">
              <div class="no-results js-not-filterable">Nothing to show</div>
                <div class="commitish-item branch-commitish selector-item js-navigation-item js-navigation-target">
                  <h4>
                      <a href="/yrosseel/lavaan/blob/master/R/02lavaanUser.R" class="js-navigation-open" data-name="master" rel="nofollow">master</a>
                  </h4>
                </div>
                <div class="commitish-item branch-commitish selector-item js-navigation-item js-navigation-target">
                  <h4>
                      <a href="/yrosseel/lavaan/blob/official/R/02lavaanUser.R" class="js-navigation-open" data-name="official" rel="nofollow">official</a>
                  </h4>
                </div>
            </div>

            <div class="js-filter-tab js-filter-tags" style="display:none" data-filterable-for="context-commitish-filter-field">
              <div class="no-results js-not-filterable">Nothing to show</div>
            </div>
          </div>
        </div><!-- /.commitish-context-context -->
      </div>

    </li>
  </ul>

  <ul class="subnav with-scope">

    <li><a href="/yrosseel/lavaan" class="selected" highlight="repo_source">Files</a></li>
    <li><a href="/yrosseel/lavaan/commits/master" highlight="repo_commits">Commits</a></li>
    <li><a href="/yrosseel/lavaan/branches" class="" highlight="repo_branches" rel="nofollow">Branches <span class="counter">2</span></a></li>
  </ul>

</div>

  
  
  


          

        </div><!-- /.repohead -->

        





<!-- block_view_fragment_key: views10/v8/blob:v21:ed26842bd1acb306d643dc3b7c0c660d -->
  <div id="slider">

    <div class="breadcrumb" data-path="R/02lavaanUser.R/">
      <b itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/yrosseel/lavaan/tree/79b51b6496090dca20de417c2b60c57d3f0c93d2" class="js-rewrite-sha" itemprop="url"><span itemprop="title">lavaan</span></a></b> / <span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/yrosseel/lavaan/tree/79b51b6496090dca20de417c2b60c57d3f0c93d2/R" class="js-rewrite-sha" itemscope="url"><span itemprop="title">R</span></a></span> / <strong class="final-path">02lavaanUser.R</strong> <span class="js-clippy mini-icon mini-icon-clippy " data-clipboard-text="R/02lavaanUser.R" data-copied-hint="copied!" data-copy-hint="copy to clipboard"></span>
    </div>


      <div class="commit file-history-tease" data-path="R/02lavaanUser.R/">
        <img class="main-avatar" height="24" src="https://secure.gravatar.com/avatar/1c4fb23eecd8db6e04f11e739df21505?s=140&amp;d=https://a248.e.akamai.net/assets.github.com%2Fimages%2Fgravatars%2Fgravatar-140.png" width="24" />
        <span class="author"><a href="/yrosseel">yrosseel</a></span>
        <time class="js-relative-date" datetime="2012-04-28T06:01:22-07:00" title="2012-04-28 06:01:22">April 28, 2012</time>
        <div class="commit-title">
            <a href="/yrosseel/lavaan/commit/ad47b18d8759ac4c7ed74c3fb94c9132abcb769b" class="message">changed name of </a><a href="https://github.com/user" class="user-mention">@user</a><a href="/yrosseel/lavaan/commit/ad47b18d8759ac4c7ed74c3fb94c9132abcb769b" class="message"> slot to @ParTable slot; same for all `user' arg…</a>
        </div>

        <div class="participation">
          <p class="quickstat"><a href="#blob_contributors_box" rel="facebox"><strong>1</strong> contributor</a></p>
          
        </div>
        <div id="blob_contributors_box" style="display:none">
          <h2>Users on GitHub who have contributed to this file</h2>
          <ul class="facebox-user-list">
            <li>
              <img height="24" src="https://secure.gravatar.com/avatar/1c4fb23eecd8db6e04f11e739df21505?s=140&amp;d=https://a248.e.akamai.net/assets.github.com%2Fimages%2Fgravatars%2Fgravatar-140.png" width="24" />
              <a href="/yrosseel">yrosseel</a>
            </li>
          </ul>
        </div>
      </div>

    <div class="frames">
      <div class="frame frame-center" data-path="R/02lavaanUser.R/" data-permalink-url="/yrosseel/lavaan/blob/79b51b6496090dca20de417c2b60c57d3f0c93d2/R/02lavaanUser.R" data-title="lavaan/R/02lavaanUser.R at master · yrosseel/lavaan · GitHub" data-type="blob">

        <div id="files" class="bubble">
          <div class="file">
            <div class="meta">
              <div class="info">
                <span class="icon"><b class="mini-icon mini-icon-text-file"></b></span>
                <span class="mode" title="File Mode">100644</span>
                  <span>709 lines (633 sloc)</span>
                <span>27.639 kb</span>
              </div>
              <ul class="button-group actions">
                  <li>
                    <a class="grouped-button file-edit-link minibutton bigger lighter js-rewrite-sha" href="/yrosseel/lavaan/edit/79b51b6496090dca20de417c2b60c57d3f0c93d2/R/02lavaanUser.R" data-method="post" rel="nofollow" data-hotkey="e">Edit this file</a>
                  </li>

                <li>
                  <a href="/yrosseel/lavaan/raw/master/R/02lavaanUser.R" class="minibutton btn-raw grouped-button bigger lighter" id="raw-url"><span class="icon"></span>Raw</a>
                </li>
                  <li>
                    <a href="/yrosseel/lavaan/blame/master/R/02lavaanUser.R" class="minibutton btn-blame grouped-button bigger lighter"><span class="icon"></span>Blame</a>
                  </li>
                <li>
                  <a href="/yrosseel/lavaan/commits/master/R/02lavaanUser.R" class="minibutton btn-history grouped-button bigger lighter" rel="nofollow"><span class="icon"></span>History</a>
                </li>
              </ul>
            </div>
              <div class="data type-r">
      <table cellpadding="0" cellspacing="0" class="lines">
        <tr>
          <td>
            <pre class="line_numbers"><span id="L1" rel="#L1">1</span>
<span id="L2" rel="#L2">2</span>
<span id="L3" rel="#L3">3</span>
<span id="L4" rel="#L4">4</span>
<span id="L5" rel="#L5">5</span>
<span id="L6" rel="#L6">6</span>
<span id="L7" rel="#L7">7</span>
<span id="L8" rel="#L8">8</span>
<span id="L9" rel="#L9">9</span>
<span id="L10" rel="#L10">10</span>
<span id="L11" rel="#L11">11</span>
<span id="L12" rel="#L12">12</span>
<span id="L13" rel="#L13">13</span>
<span id="L14" rel="#L14">14</span>
<span id="L15" rel="#L15">15</span>
<span id="L16" rel="#L16">16</span>
<span id="L17" rel="#L17">17</span>
<span id="L18" rel="#L18">18</span>
<span id="L19" rel="#L19">19</span>
<span id="L20" rel="#L20">20</span>
<span id="L21" rel="#L21">21</span>
<span id="L22" rel="#L22">22</span>
<span id="L23" rel="#L23">23</span>
<span id="L24" rel="#L24">24</span>
<span id="L25" rel="#L25">25</span>
<span id="L26" rel="#L26">26</span>
<span id="L27" rel="#L27">27</span>
<span id="L28" rel="#L28">28</span>
<span id="L29" rel="#L29">29</span>
<span id="L30" rel="#L30">30</span>
<span id="L31" rel="#L31">31</span>
<span id="L32" rel="#L32">32</span>
<span id="L33" rel="#L33">33</span>
<span id="L34" rel="#L34">34</span>
<span id="L35" rel="#L35">35</span>
<span id="L36" rel="#L36">36</span>
<span id="L37" rel="#L37">37</span>
<span id="L38" rel="#L38">38</span>
<span id="L39" rel="#L39">39</span>
<span id="L40" rel="#L40">40</span>
<span id="L41" rel="#L41">41</span>
<span id="L42" rel="#L42">42</span>
<span id="L43" rel="#L43">43</span>
<span id="L44" rel="#L44">44</span>
<span id="L45" rel="#L45">45</span>
<span id="L46" rel="#L46">46</span>
<span id="L47" rel="#L47">47</span>
<span id="L48" rel="#L48">48</span>
<span id="L49" rel="#L49">49</span>
<span id="L50" rel="#L50">50</span>
<span id="L51" rel="#L51">51</span>
<span id="L52" rel="#L52">52</span>
<span id="L53" rel="#L53">53</span>
<span id="L54" rel="#L54">54</span>
<span id="L55" rel="#L55">55</span>
<span id="L56" rel="#L56">56</span>
<span id="L57" rel="#L57">57</span>
<span id="L58" rel="#L58">58</span>
<span id="L59" rel="#L59">59</span>
<span id="L60" rel="#L60">60</span>
<span id="L61" rel="#L61">61</span>
<span id="L62" rel="#L62">62</span>
<span id="L63" rel="#L63">63</span>
<span id="L64" rel="#L64">64</span>
<span id="L65" rel="#L65">65</span>
<span id="L66" rel="#L66">66</span>
<span id="L67" rel="#L67">67</span>
<span id="L68" rel="#L68">68</span>
<span id="L69" rel="#L69">69</span>
<span id="L70" rel="#L70">70</span>
<span id="L71" rel="#L71">71</span>
<span id="L72" rel="#L72">72</span>
<span id="L73" rel="#L73">73</span>
<span id="L74" rel="#L74">74</span>
<span id="L75" rel="#L75">75</span>
<span id="L76" rel="#L76">76</span>
<span id="L77" rel="#L77">77</span>
<span id="L78" rel="#L78">78</span>
<span id="L79" rel="#L79">79</span>
<span id="L80" rel="#L80">80</span>
<span id="L81" rel="#L81">81</span>
<span id="L82" rel="#L82">82</span>
<span id="L83" rel="#L83">83</span>
<span id="L84" rel="#L84">84</span>
<span id="L85" rel="#L85">85</span>
<span id="L86" rel="#L86">86</span>
<span id="L87" rel="#L87">87</span>
<span id="L88" rel="#L88">88</span>
<span id="L89" rel="#L89">89</span>
<span id="L90" rel="#L90">90</span>
<span id="L91" rel="#L91">91</span>
<span id="L92" rel="#L92">92</span>
<span id="L93" rel="#L93">93</span>
<span id="L94" rel="#L94">94</span>
<span id="L95" rel="#L95">95</span>
<span id="L96" rel="#L96">96</span>
<span id="L97" rel="#L97">97</span>
<span id="L98" rel="#L98">98</span>
<span id="L99" rel="#L99">99</span>
<span id="L100" rel="#L100">100</span>
<span id="L101" rel="#L101">101</span>
<span id="L102" rel="#L102">102</span>
<span id="L103" rel="#L103">103</span>
<span id="L104" rel="#L104">104</span>
<span id="L105" rel="#L105">105</span>
<span id="L106" rel="#L106">106</span>
<span id="L107" rel="#L107">107</span>
<span id="L108" rel="#L108">108</span>
<span id="L109" rel="#L109">109</span>
<span id="L110" rel="#L110">110</span>
<span id="L111" rel="#L111">111</span>
<span id="L112" rel="#L112">112</span>
<span id="L113" rel="#L113">113</span>
<span id="L114" rel="#L114">114</span>
<span id="L115" rel="#L115">115</span>
<span id="L116" rel="#L116">116</span>
<span id="L117" rel="#L117">117</span>
<span id="L118" rel="#L118">118</span>
<span id="L119" rel="#L119">119</span>
<span id="L120" rel="#L120">120</span>
<span id="L121" rel="#L121">121</span>
<span id="L122" rel="#L122">122</span>
<span id="L123" rel="#L123">123</span>
<span id="L124" rel="#L124">124</span>
<span id="L125" rel="#L125">125</span>
<span id="L126" rel="#L126">126</span>
<span id="L127" rel="#L127">127</span>
<span id="L128" rel="#L128">128</span>
<span id="L129" rel="#L129">129</span>
<span id="L130" rel="#L130">130</span>
<span id="L131" rel="#L131">131</span>
<span id="L132" rel="#L132">132</span>
<span id="L133" rel="#L133">133</span>
<span id="L134" rel="#L134">134</span>
<span id="L135" rel="#L135">135</span>
<span id="L136" rel="#L136">136</span>
<span id="L137" rel="#L137">137</span>
<span id="L138" rel="#L138">138</span>
<span id="L139" rel="#L139">139</span>
<span id="L140" rel="#L140">140</span>
<span id="L141" rel="#L141">141</span>
<span id="L142" rel="#L142">142</span>
<span id="L143" rel="#L143">143</span>
<span id="L144" rel="#L144">144</span>
<span id="L145" rel="#L145">145</span>
<span id="L146" rel="#L146">146</span>
<span id="L147" rel="#L147">147</span>
<span id="L148" rel="#L148">148</span>
<span id="L149" rel="#L149">149</span>
<span id="L150" rel="#L150">150</span>
<span id="L151" rel="#L151">151</span>
<span id="L152" rel="#L152">152</span>
<span id="L153" rel="#L153">153</span>
<span id="L154" rel="#L154">154</span>
<span id="L155" rel="#L155">155</span>
<span id="L156" rel="#L156">156</span>
<span id="L157" rel="#L157">157</span>
<span id="L158" rel="#L158">158</span>
<span id="L159" rel="#L159">159</span>
<span id="L160" rel="#L160">160</span>
<span id="L161" rel="#L161">161</span>
<span id="L162" rel="#L162">162</span>
<span id="L163" rel="#L163">163</span>
<span id="L164" rel="#L164">164</span>
<span id="L165" rel="#L165">165</span>
<span id="L166" rel="#L166">166</span>
<span id="L167" rel="#L167">167</span>
<span id="L168" rel="#L168">168</span>
<span id="L169" rel="#L169">169</span>
<span id="L170" rel="#L170">170</span>
<span id="L171" rel="#L171">171</span>
<span id="L172" rel="#L172">172</span>
<span id="L173" rel="#L173">173</span>
<span id="L174" rel="#L174">174</span>
<span id="L175" rel="#L175">175</span>
<span id="L176" rel="#L176">176</span>
<span id="L177" rel="#L177">177</span>
<span id="L178" rel="#L178">178</span>
<span id="L179" rel="#L179">179</span>
<span id="L180" rel="#L180">180</span>
<span id="L181" rel="#L181">181</span>
<span id="L182" rel="#L182">182</span>
<span id="L183" rel="#L183">183</span>
<span id="L184" rel="#L184">184</span>
<span id="L185" rel="#L185">185</span>
<span id="L186" rel="#L186">186</span>
<span id="L187" rel="#L187">187</span>
<span id="L188" rel="#L188">188</span>
<span id="L189" rel="#L189">189</span>
<span id="L190" rel="#L190">190</span>
<span id="L191" rel="#L191">191</span>
<span id="L192" rel="#L192">192</span>
<span id="L193" rel="#L193">193</span>
<span id="L194" rel="#L194">194</span>
<span id="L195" rel="#L195">195</span>
<span id="L196" rel="#L196">196</span>
<span id="L197" rel="#L197">197</span>
<span id="L198" rel="#L198">198</span>
<span id="L199" rel="#L199">199</span>
<span id="L200" rel="#L200">200</span>
<span id="L201" rel="#L201">201</span>
<span id="L202" rel="#L202">202</span>
<span id="L203" rel="#L203">203</span>
<span id="L204" rel="#L204">204</span>
<span id="L205" rel="#L205">205</span>
<span id="L206" rel="#L206">206</span>
<span id="L207" rel="#L207">207</span>
<span id="L208" rel="#L208">208</span>
<span id="L209" rel="#L209">209</span>
<span id="L210" rel="#L210">210</span>
<span id="L211" rel="#L211">211</span>
<span id="L212" rel="#L212">212</span>
<span id="L213" rel="#L213">213</span>
<span id="L214" rel="#L214">214</span>
<span id="L215" rel="#L215">215</span>
<span id="L216" rel="#L216">216</span>
<span id="L217" rel="#L217">217</span>
<span id="L218" rel="#L218">218</span>
<span id="L219" rel="#L219">219</span>
<span id="L220" rel="#L220">220</span>
<span id="L221" rel="#L221">221</span>
<span id="L222" rel="#L222">222</span>
<span id="L223" rel="#L223">223</span>
<span id="L224" rel="#L224">224</span>
<span id="L225" rel="#L225">225</span>
<span id="L226" rel="#L226">226</span>
<span id="L227" rel="#L227">227</span>
<span id="L228" rel="#L228">228</span>
<span id="L229" rel="#L229">229</span>
<span id="L230" rel="#L230">230</span>
<span id="L231" rel="#L231">231</span>
<span id="L232" rel="#L232">232</span>
<span id="L233" rel="#L233">233</span>
<span id="L234" rel="#L234">234</span>
<span id="L235" rel="#L235">235</span>
<span id="L236" rel="#L236">236</span>
<span id="L237" rel="#L237">237</span>
<span id="L238" rel="#L238">238</span>
<span id="L239" rel="#L239">239</span>
<span id="L240" rel="#L240">240</span>
<span id="L241" rel="#L241">241</span>
<span id="L242" rel="#L242">242</span>
<span id="L243" rel="#L243">243</span>
<span id="L244" rel="#L244">244</span>
<span id="L245" rel="#L245">245</span>
<span id="L246" rel="#L246">246</span>
<span id="L247" rel="#L247">247</span>
<span id="L248" rel="#L248">248</span>
<span id="L249" rel="#L249">249</span>
<span id="L250" rel="#L250">250</span>
<span id="L251" rel="#L251">251</span>
<span id="L252" rel="#L252">252</span>
<span id="L253" rel="#L253">253</span>
<span id="L254" rel="#L254">254</span>
<span id="L255" rel="#L255">255</span>
<span id="L256" rel="#L256">256</span>
<span id="L257" rel="#L257">257</span>
<span id="L258" rel="#L258">258</span>
<span id="L259" rel="#L259">259</span>
<span id="L260" rel="#L260">260</span>
<span id="L261" rel="#L261">261</span>
<span id="L262" rel="#L262">262</span>
<span id="L263" rel="#L263">263</span>
<span id="L264" rel="#L264">264</span>
<span id="L265" rel="#L265">265</span>
<span id="L266" rel="#L266">266</span>
<span id="L267" rel="#L267">267</span>
<span id="L268" rel="#L268">268</span>
<span id="L269" rel="#L269">269</span>
<span id="L270" rel="#L270">270</span>
<span id="L271" rel="#L271">271</span>
<span id="L272" rel="#L272">272</span>
<span id="L273" rel="#L273">273</span>
<span id="L274" rel="#L274">274</span>
<span id="L275" rel="#L275">275</span>
<span id="L276" rel="#L276">276</span>
<span id="L277" rel="#L277">277</span>
<span id="L278" rel="#L278">278</span>
<span id="L279" rel="#L279">279</span>
<span id="L280" rel="#L280">280</span>
<span id="L281" rel="#L281">281</span>
<span id="L282" rel="#L282">282</span>
<span id="L283" rel="#L283">283</span>
<span id="L284" rel="#L284">284</span>
<span id="L285" rel="#L285">285</span>
<span id="L286" rel="#L286">286</span>
<span id="L287" rel="#L287">287</span>
<span id="L288" rel="#L288">288</span>
<span id="L289" rel="#L289">289</span>
<span id="L290" rel="#L290">290</span>
<span id="L291" rel="#L291">291</span>
<span id="L292" rel="#L292">292</span>
<span id="L293" rel="#L293">293</span>
<span id="L294" rel="#L294">294</span>
<span id="L295" rel="#L295">295</span>
<span id="L296" rel="#L296">296</span>
<span id="L297" rel="#L297">297</span>
<span id="L298" rel="#L298">298</span>
<span id="L299" rel="#L299">299</span>
<span id="L300" rel="#L300">300</span>
<span id="L301" rel="#L301">301</span>
<span id="L302" rel="#L302">302</span>
<span id="L303" rel="#L303">303</span>
<span id="L304" rel="#L304">304</span>
<span id="L305" rel="#L305">305</span>
<span id="L306" rel="#L306">306</span>
<span id="L307" rel="#L307">307</span>
<span id="L308" rel="#L308">308</span>
<span id="L309" rel="#L309">309</span>
<span id="L310" rel="#L310">310</span>
<span id="L311" rel="#L311">311</span>
<span id="L312" rel="#L312">312</span>
<span id="L313" rel="#L313">313</span>
<span id="L314" rel="#L314">314</span>
<span id="L315" rel="#L315">315</span>
<span id="L316" rel="#L316">316</span>
<span id="L317" rel="#L317">317</span>
<span id="L318" rel="#L318">318</span>
<span id="L319" rel="#L319">319</span>
<span id="L320" rel="#L320">320</span>
<span id="L321" rel="#L321">321</span>
<span id="L322" rel="#L322">322</span>
<span id="L323" rel="#L323">323</span>
<span id="L324" rel="#L324">324</span>
<span id="L325" rel="#L325">325</span>
<span id="L326" rel="#L326">326</span>
<span id="L327" rel="#L327">327</span>
<span id="L328" rel="#L328">328</span>
<span id="L329" rel="#L329">329</span>
<span id="L330" rel="#L330">330</span>
<span id="L331" rel="#L331">331</span>
<span id="L332" rel="#L332">332</span>
<span id="L333" rel="#L333">333</span>
<span id="L334" rel="#L334">334</span>
<span id="L335" rel="#L335">335</span>
<span id="L336" rel="#L336">336</span>
<span id="L337" rel="#L337">337</span>
<span id="L338" rel="#L338">338</span>
<span id="L339" rel="#L339">339</span>
<span id="L340" rel="#L340">340</span>
<span id="L341" rel="#L341">341</span>
<span id="L342" rel="#L342">342</span>
<span id="L343" rel="#L343">343</span>
<span id="L344" rel="#L344">344</span>
<span id="L345" rel="#L345">345</span>
<span id="L346" rel="#L346">346</span>
<span id="L347" rel="#L347">347</span>
<span id="L348" rel="#L348">348</span>
<span id="L349" rel="#L349">349</span>
<span id="L350" rel="#L350">350</span>
<span id="L351" rel="#L351">351</span>
<span id="L352" rel="#L352">352</span>
<span id="L353" rel="#L353">353</span>
<span id="L354" rel="#L354">354</span>
<span id="L355" rel="#L355">355</span>
<span id="L356" rel="#L356">356</span>
<span id="L357" rel="#L357">357</span>
<span id="L358" rel="#L358">358</span>
<span id="L359" rel="#L359">359</span>
<span id="L360" rel="#L360">360</span>
<span id="L361" rel="#L361">361</span>
<span id="L362" rel="#L362">362</span>
<span id="L363" rel="#L363">363</span>
<span id="L364" rel="#L364">364</span>
<span id="L365" rel="#L365">365</span>
<span id="L366" rel="#L366">366</span>
<span id="L367" rel="#L367">367</span>
<span id="L368" rel="#L368">368</span>
<span id="L369" rel="#L369">369</span>
<span id="L370" rel="#L370">370</span>
<span id="L371" rel="#L371">371</span>
<span id="L372" rel="#L372">372</span>
<span id="L373" rel="#L373">373</span>
<span id="L374" rel="#L374">374</span>
<span id="L375" rel="#L375">375</span>
<span id="L376" rel="#L376">376</span>
<span id="L377" rel="#L377">377</span>
<span id="L378" rel="#L378">378</span>
<span id="L379" rel="#L379">379</span>
<span id="L380" rel="#L380">380</span>
<span id="L381" rel="#L381">381</span>
<span id="L382" rel="#L382">382</span>
<span id="L383" rel="#L383">383</span>
<span id="L384" rel="#L384">384</span>
<span id="L385" rel="#L385">385</span>
<span id="L386" rel="#L386">386</span>
<span id="L387" rel="#L387">387</span>
<span id="L388" rel="#L388">388</span>
<span id="L389" rel="#L389">389</span>
<span id="L390" rel="#L390">390</span>
<span id="L391" rel="#L391">391</span>
<span id="L392" rel="#L392">392</span>
<span id="L393" rel="#L393">393</span>
<span id="L394" rel="#L394">394</span>
<span id="L395" rel="#L395">395</span>
<span id="L396" rel="#L396">396</span>
<span id="L397" rel="#L397">397</span>
<span id="L398" rel="#L398">398</span>
<span id="L399" rel="#L399">399</span>
<span id="L400" rel="#L400">400</span>
<span id="L401" rel="#L401">401</span>
<span id="L402" rel="#L402">402</span>
<span id="L403" rel="#L403">403</span>
<span id="L404" rel="#L404">404</span>
<span id="L405" rel="#L405">405</span>
<span id="L406" rel="#L406">406</span>
<span id="L407" rel="#L407">407</span>
<span id="L408" rel="#L408">408</span>
<span id="L409" rel="#L409">409</span>
<span id="L410" rel="#L410">410</span>
<span id="L411" rel="#L411">411</span>
<span id="L412" rel="#L412">412</span>
<span id="L413" rel="#L413">413</span>
<span id="L414" rel="#L414">414</span>
<span id="L415" rel="#L415">415</span>
<span id="L416" rel="#L416">416</span>
<span id="L417" rel="#L417">417</span>
<span id="L418" rel="#L418">418</span>
<span id="L419" rel="#L419">419</span>
<span id="L420" rel="#L420">420</span>
<span id="L421" rel="#L421">421</span>
<span id="L422" rel="#L422">422</span>
<span id="L423" rel="#L423">423</span>
<span id="L424" rel="#L424">424</span>
<span id="L425" rel="#L425">425</span>
<span id="L426" rel="#L426">426</span>
<span id="L427" rel="#L427">427</span>
<span id="L428" rel="#L428">428</span>
<span id="L429" rel="#L429">429</span>
<span id="L430" rel="#L430">430</span>
<span id="L431" rel="#L431">431</span>
<span id="L432" rel="#L432">432</span>
<span id="L433" rel="#L433">433</span>
<span id="L434" rel="#L434">434</span>
<span id="L435" rel="#L435">435</span>
<span id="L436" rel="#L436">436</span>
<span id="L437" rel="#L437">437</span>
<span id="L438" rel="#L438">438</span>
<span id="L439" rel="#L439">439</span>
<span id="L440" rel="#L440">440</span>
<span id="L441" rel="#L441">441</span>
<span id="L442" rel="#L442">442</span>
<span id="L443" rel="#L443">443</span>
<span id="L444" rel="#L444">444</span>
<span id="L445" rel="#L445">445</span>
<span id="L446" rel="#L446">446</span>
<span id="L447" rel="#L447">447</span>
<span id="L448" rel="#L448">448</span>
<span id="L449" rel="#L449">449</span>
<span id="L450" rel="#L450">450</span>
<span id="L451" rel="#L451">451</span>
<span id="L452" rel="#L452">452</span>
<span id="L453" rel="#L453">453</span>
<span id="L454" rel="#L454">454</span>
<span id="L455" rel="#L455">455</span>
<span id="L456" rel="#L456">456</span>
<span id="L457" rel="#L457">457</span>
<span id="L458" rel="#L458">458</span>
<span id="L459" rel="#L459">459</span>
<span id="L460" rel="#L460">460</span>
<span id="L461" rel="#L461">461</span>
<span id="L462" rel="#L462">462</span>
<span id="L463" rel="#L463">463</span>
<span id="L464" rel="#L464">464</span>
<span id="L465" rel="#L465">465</span>
<span id="L466" rel="#L466">466</span>
<span id="L467" rel="#L467">467</span>
<span id="L468" rel="#L468">468</span>
<span id="L469" rel="#L469">469</span>
<span id="L470" rel="#L470">470</span>
<span id="L471" rel="#L471">471</span>
<span id="L472" rel="#L472">472</span>
<span id="L473" rel="#L473">473</span>
<span id="L474" rel="#L474">474</span>
<span id="L475" rel="#L475">475</span>
<span id="L476" rel="#L476">476</span>
<span id="L477" rel="#L477">477</span>
<span id="L478" rel="#L478">478</span>
<span id="L479" rel="#L479">479</span>
<span id="L480" rel="#L480">480</span>
<span id="L481" rel="#L481">481</span>
<span id="L482" rel="#L482">482</span>
<span id="L483" rel="#L483">483</span>
<span id="L484" rel="#L484">484</span>
<span id="L485" rel="#L485">485</span>
<span id="L486" rel="#L486">486</span>
<span id="L487" rel="#L487">487</span>
<span id="L488" rel="#L488">488</span>
<span id="L489" rel="#L489">489</span>
<span id="L490" rel="#L490">490</span>
<span id="L491" rel="#L491">491</span>
<span id="L492" rel="#L492">492</span>
<span id="L493" rel="#L493">493</span>
<span id="L494" rel="#L494">494</span>
<span id="L495" rel="#L495">495</span>
<span id="L496" rel="#L496">496</span>
<span id="L497" rel="#L497">497</span>
<span id="L498" rel="#L498">498</span>
<span id="L499" rel="#L499">499</span>
<span id="L500" rel="#L500">500</span>
<span id="L501" rel="#L501">501</span>
<span id="L502" rel="#L502">502</span>
<span id="L503" rel="#L503">503</span>
<span id="L504" rel="#L504">504</span>
<span id="L505" rel="#L505">505</span>
<span id="L506" rel="#L506">506</span>
<span id="L507" rel="#L507">507</span>
<span id="L508" rel="#L508">508</span>
<span id="L509" rel="#L509">509</span>
<span id="L510" rel="#L510">510</span>
<span id="L511" rel="#L511">511</span>
<span id="L512" rel="#L512">512</span>
<span id="L513" rel="#L513">513</span>
<span id="L514" rel="#L514">514</span>
<span id="L515" rel="#L515">515</span>
<span id="L516" rel="#L516">516</span>
<span id="L517" rel="#L517">517</span>
<span id="L518" rel="#L518">518</span>
<span id="L519" rel="#L519">519</span>
<span id="L520" rel="#L520">520</span>
<span id="L521" rel="#L521">521</span>
<span id="L522" rel="#L522">522</span>
<span id="L523" rel="#L523">523</span>
<span id="L524" rel="#L524">524</span>
<span id="L525" rel="#L525">525</span>
<span id="L526" rel="#L526">526</span>
<span id="L527" rel="#L527">527</span>
<span id="L528" rel="#L528">528</span>
<span id="L529" rel="#L529">529</span>
<span id="L530" rel="#L530">530</span>
<span id="L531" rel="#L531">531</span>
<span id="L532" rel="#L532">532</span>
<span id="L533" rel="#L533">533</span>
<span id="L534" rel="#L534">534</span>
<span id="L535" rel="#L535">535</span>
<span id="L536" rel="#L536">536</span>
<span id="L537" rel="#L537">537</span>
<span id="L538" rel="#L538">538</span>
<span id="L539" rel="#L539">539</span>
<span id="L540" rel="#L540">540</span>
<span id="L541" rel="#L541">541</span>
<span id="L542" rel="#L542">542</span>
<span id="L543" rel="#L543">543</span>
<span id="L544" rel="#L544">544</span>
<span id="L545" rel="#L545">545</span>
<span id="L546" rel="#L546">546</span>
<span id="L547" rel="#L547">547</span>
<span id="L548" rel="#L548">548</span>
<span id="L549" rel="#L549">549</span>
<span id="L550" rel="#L550">550</span>
<span id="L551" rel="#L551">551</span>
<span id="L552" rel="#L552">552</span>
<span id="L553" rel="#L553">553</span>
<span id="L554" rel="#L554">554</span>
<span id="L555" rel="#L555">555</span>
<span id="L556" rel="#L556">556</span>
<span id="L557" rel="#L557">557</span>
<span id="L558" rel="#L558">558</span>
<span id="L559" rel="#L559">559</span>
<span id="L560" rel="#L560">560</span>
<span id="L561" rel="#L561">561</span>
<span id="L562" rel="#L562">562</span>
<span id="L563" rel="#L563">563</span>
<span id="L564" rel="#L564">564</span>
<span id="L565" rel="#L565">565</span>
<span id="L566" rel="#L566">566</span>
<span id="L567" rel="#L567">567</span>
<span id="L568" rel="#L568">568</span>
<span id="L569" rel="#L569">569</span>
<span id="L570" rel="#L570">570</span>
<span id="L571" rel="#L571">571</span>
<span id="L572" rel="#L572">572</span>
<span id="L573" rel="#L573">573</span>
<span id="L574" rel="#L574">574</span>
<span id="L575" rel="#L575">575</span>
<span id="L576" rel="#L576">576</span>
<span id="L577" rel="#L577">577</span>
<span id="L578" rel="#L578">578</span>
<span id="L579" rel="#L579">579</span>
<span id="L580" rel="#L580">580</span>
<span id="L581" rel="#L581">581</span>
<span id="L582" rel="#L582">582</span>
<span id="L583" rel="#L583">583</span>
<span id="L584" rel="#L584">584</span>
<span id="L585" rel="#L585">585</span>
<span id="L586" rel="#L586">586</span>
<span id="L587" rel="#L587">587</span>
<span id="L588" rel="#L588">588</span>
<span id="L589" rel="#L589">589</span>
<span id="L590" rel="#L590">590</span>
<span id="L591" rel="#L591">591</span>
<span id="L592" rel="#L592">592</span>
<span id="L593" rel="#L593">593</span>
<span id="L594" rel="#L594">594</span>
<span id="L595" rel="#L595">595</span>
<span id="L596" rel="#L596">596</span>
<span id="L597" rel="#L597">597</span>
<span id="L598" rel="#L598">598</span>
<span id="L599" rel="#L599">599</span>
<span id="L600" rel="#L600">600</span>
<span id="L601" rel="#L601">601</span>
<span id="L602" rel="#L602">602</span>
<span id="L603" rel="#L603">603</span>
<span id="L604" rel="#L604">604</span>
<span id="L605" rel="#L605">605</span>
<span id="L606" rel="#L606">606</span>
<span id="L607" rel="#L607">607</span>
<span id="L608" rel="#L608">608</span>
<span id="L609" rel="#L609">609</span>
<span id="L610" rel="#L610">610</span>
<span id="L611" rel="#L611">611</span>
<span id="L612" rel="#L612">612</span>
<span id="L613" rel="#L613">613</span>
<span id="L614" rel="#L614">614</span>
<span id="L615" rel="#L615">615</span>
<span id="L616" rel="#L616">616</span>
<span id="L617" rel="#L617">617</span>
<span id="L618" rel="#L618">618</span>
<span id="L619" rel="#L619">619</span>
<span id="L620" rel="#L620">620</span>
<span id="L621" rel="#L621">621</span>
<span id="L622" rel="#L622">622</span>
<span id="L623" rel="#L623">623</span>
<span id="L624" rel="#L624">624</span>
<span id="L625" rel="#L625">625</span>
<span id="L626" rel="#L626">626</span>
<span id="L627" rel="#L627">627</span>
<span id="L628" rel="#L628">628</span>
<span id="L629" rel="#L629">629</span>
<span id="L630" rel="#L630">630</span>
<span id="L631" rel="#L631">631</span>
<span id="L632" rel="#L632">632</span>
<span id="L633" rel="#L633">633</span>
<span id="L634" rel="#L634">634</span>
<span id="L635" rel="#L635">635</span>
<span id="L636" rel="#L636">636</span>
<span id="L637" rel="#L637">637</span>
<span id="L638" rel="#L638">638</span>
<span id="L639" rel="#L639">639</span>
<span id="L640" rel="#L640">640</span>
<span id="L641" rel="#L641">641</span>
<span id="L642" rel="#L642">642</span>
<span id="L643" rel="#L643">643</span>
<span id="L644" rel="#L644">644</span>
<span id="L645" rel="#L645">645</span>
<span id="L646" rel="#L646">646</span>
<span id="L647" rel="#L647">647</span>
<span id="L648" rel="#L648">648</span>
<span id="L649" rel="#L649">649</span>
<span id="L650" rel="#L650">650</span>
<span id="L651" rel="#L651">651</span>
<span id="L652" rel="#L652">652</span>
<span id="L653" rel="#L653">653</span>
<span id="L654" rel="#L654">654</span>
<span id="L655" rel="#L655">655</span>
<span id="L656" rel="#L656">656</span>
<span id="L657" rel="#L657">657</span>
<span id="L658" rel="#L658">658</span>
<span id="L659" rel="#L659">659</span>
<span id="L660" rel="#L660">660</span>
<span id="L661" rel="#L661">661</span>
<span id="L662" rel="#L662">662</span>
<span id="L663" rel="#L663">663</span>
<span id="L664" rel="#L664">664</span>
<span id="L665" rel="#L665">665</span>
<span id="L666" rel="#L666">666</span>
<span id="L667" rel="#L667">667</span>
<span id="L668" rel="#L668">668</span>
<span id="L669" rel="#L669">669</span>
<span id="L670" rel="#L670">670</span>
<span id="L671" rel="#L671">671</span>
<span id="L672" rel="#L672">672</span>
<span id="L673" rel="#L673">673</span>
<span id="L674" rel="#L674">674</span>
<span id="L675" rel="#L675">675</span>
<span id="L676" rel="#L676">676</span>
<span id="L677" rel="#L677">677</span>
<span id="L678" rel="#L678">678</span>
<span id="L679" rel="#L679">679</span>
<span id="L680" rel="#L680">680</span>
<span id="L681" rel="#L681">681</span>
<span id="L682" rel="#L682">682</span>
<span id="L683" rel="#L683">683</span>
<span id="L684" rel="#L684">684</span>
<span id="L685" rel="#L685">685</span>
<span id="L686" rel="#L686">686</span>
<span id="L687" rel="#L687">687</span>
<span id="L688" rel="#L688">688</span>
<span id="L689" rel="#L689">689</span>
<span id="L690" rel="#L690">690</span>
<span id="L691" rel="#L691">691</span>
<span id="L692" rel="#L692">692</span>
<span id="L693" rel="#L693">693</span>
<span id="L694" rel="#L694">694</span>
<span id="L695" rel="#L695">695</span>
<span id="L696" rel="#L696">696</span>
<span id="L697" rel="#L697">697</span>
<span id="L698" rel="#L698">698</span>
<span id="L699" rel="#L699">699</span>
<span id="L700" rel="#L700">700</span>
<span id="L701" rel="#L701">701</span>
<span id="L702" rel="#L702">702</span>
<span id="L703" rel="#L703">703</span>
<span id="L704" rel="#L704">704</span>
<span id="L705" rel="#L705">705</span>
<span id="L706" rel="#L706">706</span>
<span id="L707" rel="#L707">707</span>
<span id="L708" rel="#L708">708</span>
</pre>
          </td>
          <td width="100%">
                <div class="highlight"><pre><div class='line' id='LC1'><span class="c1"># constructor for the lavParTable model description</span></div><div class='line' id='LC2'><span class="c1">#</span></div><div class='line' id='LC3'><span class="c1"># initial version: YR 22/05/2009</span></div><div class='line' id='LC4'><span class="c1">#  major revision: YR 02/11/2010: - FLATTEN the model syntax and turn it into a</span></div><div class='line' id='LC5'><span class="c1">#                                   data.frame, with a &quot;modifiers&quot; attribute</span></div><div class='line' id='LC6'><span class="c1">#                                 - add default elements here</span></div><div class='line' id='LC7'><span class="c1">#                                 - check for duplicate elements</span></div><div class='line' id='LC8'><span class="c1">#                                 - allow for every possible model...</span></div><div class='line' id='LC9'><span class="c1">#                                 - since 0.4-5</span></div><div class='line' id='LC10'><span class="c1">#                                 - the end result is a full description of</span></div><div class='line' id='LC11'><span class="c1">#                                   a model (but no matrix representation)</span></div><div class='line' id='LC12'><br/></div><div class='line' id='LC13'><br/></div><div class='line' id='LC14'>lavaanify <span class="o">&lt;-</span> lavParTable <span class="o">&lt;-</span> <span class="kr">function</span><span class="p">(</span></div><div class='line' id='LC15'><br/></div><div class='line' id='LC16'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;model           <span class="o">=</span> <span class="kc">NULL</span><span class="p">,</span> </div><div class='line' id='LC17'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;meanstructure   <span class="o">=</span> <span class="kc">FALSE</span><span class="p">,</span></div><div class='line' id='LC18'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;int.ov.free     <span class="o">=</span> <span class="kc">FALSE</span><span class="p">,</span></div><div class='line' id='LC19'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;int.lv.free     <span class="o">=</span> <span class="kc">FALSE</span><span class="p">,</span></div><div class='line' id='LC20'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;orthogonal      <span class="o">=</span> <span class="kc">FALSE</span><span class="p">,</span> </div><div class='line' id='LC21'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;std.lv          <span class="o">=</span> <span class="kc">FALSE</span><span class="p">,</span></div><div class='line' id='LC22'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;fixed.x         <span class="o">=</span> <span class="kc">TRUE</span><span class="p">,</span></div><div class='line' id='LC23'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;constraints     <span class="o">=</span> <span class="kc">NULL</span><span class="p">,</span></div><div class='line' id='LC24'><br/></div><div class='line' id='LC25'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;auto            <span class="o">=</span> <span class="kc">FALSE</span><span class="p">,</span></div><div class='line' id='LC26'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;model.type      <span class="o">=</span> <span class="s">&quot;sem&quot;</span><span class="p">,</span></div><div class='line' id='LC27'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;auto.fix.first  <span class="o">=</span> <span class="kc">FALSE</span><span class="p">,</span></div><div class='line' id='LC28'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;auto.fix.single <span class="o">=</span> <span class="kc">FALSE</span><span class="p">,</span></div><div class='line' id='LC29'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;auto.var        <span class="o">=</span> <span class="kc">FALSE</span><span class="p">,</span></div><div class='line' id='LC30'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;auto.cov.lv.x   <span class="o">=</span> <span class="kc">FALSE</span><span class="p">,</span></div><div class='line' id='LC31'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;auto.cov.y      <span class="o">=</span> <span class="kc">FALSE</span><span class="p">,</span></div><div class='line' id='LC32'><br/></div><div class='line' id='LC33'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;ngroups         <span class="o">=</span> <span class="m">1</span>L<span class="p">,</span></div><div class='line' id='LC34'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;group.equal     <span class="o">=</span> <span class="kc">NULL</span><span class="p">,</span></div><div class='line' id='LC35'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;group.partial   <span class="o">=</span> <span class="kc">NULL</span><span class="p">,</span></div><div class='line' id='LC36'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;debug           <span class="o">=</span> <span class="kc">FALSE</span><span class="p">,</span></div><div class='line' id='LC37'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;warn            <span class="o">=</span> <span class="kc">TRUE</span><span class="p">,</span></div><div class='line' id='LC38'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC39'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;as.data.frame.   <span class="o">=</span> <span class="kc">TRUE</span><span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC40'><br/></div><div class='line' id='LC41'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># check if model is already FLAT</span></div><div class='line' id='LC42'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>is.list<span class="p">(</span>model<span class="p">)</span> <span class="o">&amp;&amp;</span> !is.null<span class="p">(</span>model<span class="p">$</span>lhs<span class="p">))</span> <span class="p">{</span></div><div class='line' id='LC43'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;FLAT <span class="o">&lt;-</span> model</div><div class='line' id='LC44'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span> <span class="kr">else</span> <span class="p">{</span></div><div class='line' id='LC45'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># parse the model syntax and flatten the user-specified model</span></div><div class='line' id='LC46'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># return a data.frame, where each line is a model element (rhs, op, lhs)</span></div><div class='line' id='LC47'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;FLAT <span class="o">&lt;-</span> parseModelString<span class="p">(</span>model.syntax<span class="o">=</span>model<span class="p">,</span> warn<span class="o">=</span>warn<span class="p">,</span> debug<span class="o">=</span><span class="kc">FALSE</span><span class="p">)</span></div><div class='line' id='LC48'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC49'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># user-specified *modifiers* are returned as an attribute</span></div><div class='line' id='LC50'>&nbsp;&nbsp;&nbsp;&nbsp;MOD  <span class="o">&lt;-</span> attr<span class="p">(</span>FLAT<span class="p">,</span> <span class="s">&quot;modifiers&quot;</span><span class="p">);</span> attr<span class="p">(</span>FLAT<span class="p">,</span> <span class="s">&quot;modifiers&quot;</span><span class="p">)</span> <span class="o">&lt;-</span> <span class="kc">NULL</span></div><div class='line' id='LC51'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># user-specified *constraints* are returned as an attribute</span></div><div class='line' id='LC52'>&nbsp;&nbsp;&nbsp;&nbsp;CON  <span class="o">&lt;-</span> attr<span class="p">(</span>FLAT<span class="p">,</span> <span class="s">&quot;constraints&quot;</span><span class="p">);</span> attr<span class="p">(</span>FLAT<span class="p">,</span> <span class="s">&quot;constraints&quot;</span><span class="p">)</span> <span class="o">&lt;-</span> <span class="kc">NULL</span></div><div class='line' id='LC53'><br/></div><div class='line' id='LC54'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># extra constraints?</span></div><div class='line' id='LC55'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>!is.null<span class="p">(</span>constraints<span class="p">)</span> <span class="o">&amp;&amp;</span> nchar<span class="p">(</span>constraints<span class="p">)</span> <span class="o">&gt;</span> <span class="m">0</span>L<span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC56'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;FLAT2 <span class="o">&lt;-</span> parseModelString<span class="p">(</span>model.syntax<span class="o">=</span>constraints<span class="p">,</span> warn<span class="o">=</span>warn<span class="p">)</span></div><div class='line' id='LC57'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;CON2 <span class="o">&lt;-</span> attr<span class="p">(</span>FLAT2<span class="p">,</span> <span class="s">&quot;constraints&quot;</span><span class="p">);</span> rm<span class="p">(</span>FLAT2<span class="p">)</span></div><div class='line' id='LC58'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;CON <span class="o">&lt;-</span> c<span class="p">(</span>CON<span class="p">,</span> CON2<span class="p">)</span></div><div class='line' id='LC59'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC60'>&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC61'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>debug<span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC62'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;cat<span class="p">(</span><span class="s">&quot;[lavaan DEBUG]: FLAT (flattened user model):\n&quot;</span><span class="p">)</span></div><div class='line' id='LC63'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;print<span class="p">(</span>FLAT<span class="p">)</span></div><div class='line' id='LC64'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;cat<span class="p">(</span><span class="s">&quot;[lavaan DEBUG]: MOD (modifiers):\n&quot;</span><span class="p">)</span></div><div class='line' id='LC65'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;print<span class="p">(</span> str<span class="p">(</span>MOD<span class="p">)</span> <span class="p">)</span></div><div class='line' id='LC66'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;cat<span class="p">(</span><span class="s">&quot;[lavaan DEBUG]: CON (constraints):\n&quot;</span><span class="p">)</span></div><div class='line' id='LC67'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;print<span class="p">(</span> str<span class="p">(</span>CON<span class="p">)</span> <span class="p">)</span></div><div class='line' id='LC68'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC69'><br/></div><div class='line' id='LC70'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># check for wrongly specified variances/covariances/intercepts</span></div><div class='line' id='LC71'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># of exogenous variables in model syntax (if fixed.x=TRUE)</span></div><div class='line' id='LC72'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>fixed.x<span class="p">)</span> <span class="p">{</span> <span class="c1"># we ignore the groups here!</span></div><div class='line' id='LC73'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># we only call this function for the warning message</span></div><div class='line' id='LC74'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;tmp <span class="o">&lt;-</span> vnames<span class="p">(</span>FLAT<span class="p">,</span> <span class="s">&quot;ov.x&quot;</span><span class="p">,</span> warn<span class="o">=</span><span class="kc">TRUE</span><span class="p">);</span> rm<span class="p">(</span>tmp<span class="p">)</span></div><div class='line' id='LC75'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC76'><br/></div><div class='line' id='LC77'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># auto=TRUE?</span></div><div class='line' id='LC78'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>auto <span class="o">&amp;&amp;</span> model.type <span class="o">==</span> <span class="s">&quot;sem&quot;</span><span class="p">)</span> <span class="p">{</span> <span class="c1"># mimic sem/cfa auto behavior</span></div><div class='line' id='LC79'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>model.type <span class="o">==</span> <span class="s">&quot;sem&quot;</span><span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC80'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;int.ov.free     <span class="o">=</span> <span class="kc">TRUE</span></div><div class='line' id='LC81'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;int.lv.free     <span class="o">=</span> <span class="kc">FALSE</span></div><div class='line' id='LC82'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;auto.fix.first  <span class="o">=</span> !std.lv</div><div class='line' id='LC83'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;auto.fix.single <span class="o">=</span> <span class="kc">TRUE</span></div><div class='line' id='LC84'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;auto.var        <span class="o">=</span> <span class="kc">TRUE</span></div><div class='line' id='LC85'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;auto.cov.lv.x   <span class="o">=</span> <span class="kc">TRUE</span></div><div class='line' id='LC86'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;auto.cov.y      <span class="o">=</span> <span class="kc">TRUE</span></div><div class='line' id='LC87'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span> <span class="kr">else</span> </div><div class='line' id='LC88'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC89'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>model.type <span class="o">==</span> <span class="s">&quot;growth&quot;</span><span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC90'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;model.type      <span class="o">=</span> <span class="s">&quot;growth&quot;</span></div><div class='line' id='LC91'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;int.ov.free     <span class="o">=</span> <span class="kc">FALSE</span></div><div class='line' id='LC92'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;int.lv.free     <span class="o">=</span> <span class="kc">TRUE</span></div><div class='line' id='LC93'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;auto.fix.first  <span class="o">=</span> !std.lv</div><div class='line' id='LC94'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;auto.fix.single <span class="o">=</span> <span class="kc">TRUE</span></div><div class='line' id='LC95'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;auto.var        <span class="o">=</span> <span class="kc">TRUE</span></div><div class='line' id='LC96'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;auto.cov.lv.x   <span class="o">=</span> <span class="kc">TRUE</span></div><div class='line' id='LC97'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;auto.cov.y      <span class="o">=</span> <span class="kc">TRUE</span>    </div><div class='line' id='LC98'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC99'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC100'><br/></div><div class='line' id='LC101'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># check for `empty&#39; model</span></div><div class='line' id='LC102'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">#if(!any(c(FLAT$op == &quot;=~&quot;, FLAT$op == &quot;~&quot;))) {</span></div><div class='line' id='LC103'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">#    # empty model; abort or warning??</span></div><div class='line' id='LC104'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">#    #stop(&quot;lavaan ERROR: model does not contain a measurement model&quot;,</span></div><div class='line' id='LC105'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">#    #     &quot;or at least one regression formula&quot;)</span></div><div class='line' id='LC106'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">#    warning(&quot;lavaan WARNING: model does not contain a measurement model&quot;,</span></div><div class='line' id='LC107'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">#          &quot; or at least one regression formula&quot;)</span></div><div class='line' id='LC108'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">#}</span></div><div class='line' id='LC109'><br/></div><div class='line' id='LC110'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># check for meanstructure</span></div><div class='line' id='LC111'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>any<span class="p">(</span>FLAT<span class="p">$</span>op <span class="o">==</span> <span class="s">&quot;~1&quot;</span><span class="p">))</span> meanstructure <span class="o">&lt;-</span> <span class="kc">TRUE</span></div><div class='line' id='LC112'><br/></div><div class='line' id='LC113'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># check for group identifiers in the syntax (op = &quot;:&quot;)</span></div><div class='line' id='LC114'>&nbsp;&nbsp;&nbsp;&nbsp;n.group.flat <span class="o">&lt;-</span> length<span class="p">(</span>which<span class="p">(</span>FLAT<span class="p">$</span>op <span class="o">==</span> <span class="s">&quot;:&quot;</span><span class="p">))</span></div><div class='line' id='LC115'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>n.group.flat <span class="o">&gt;</span> <span class="m">0</span>L <span class="o">&amp;&amp;</span> n.group.flat <span class="o">!=</span> ngroups<span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC116'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;stop<span class="p">(</span><span class="s">&quot;lavaan ERROR: syntax defines &quot;</span><span class="p">,</span> n.group.flat<span class="p">,</span> <span class="s">&quot; groups; &quot;</span><span class="p">,</span></div><div class='line' id='LC117'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&quot;data suggests &quot;</span><span class="p">,</span> ngroups<span class="p">,</span> <span class="s">&quot; groups&quot;</span><span class="p">)</span></div><div class='line' id='LC118'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC119'><br/></div><div class='line' id='LC120'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># for each `group&#39; in n.group.flat, produce a USER</span></div><div class='line' id='LC121'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>n.group.flat <span class="o">&gt;</span> <span class="m">0</span>L<span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC122'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># split the FLAT data.frame per `group&#39;, create LIST</span></div><div class='line' id='LC123'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># for each `group&#39;, and bind them together</span></div><div class='line' id='LC124'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;group.op.idx <span class="o">&lt;-</span> c<span class="p">(</span>which<span class="p">(</span>FLAT<span class="p">$</span>op <span class="o">==</span> <span class="s">&quot;:&quot;</span><span class="p">),</span> nrow<span class="p">(</span>FLAT<span class="p">)</span><span class="o">+</span><span class="m">1</span>L<span class="p">)</span></div><div class='line' id='LC125'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">for</span><span class="p">(</span>g in <span class="m">1</span>:n.group.flat<span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC126'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;FLAT.group <span class="o">&lt;-</span> FLAT<span class="p">[(</span>group.op.idx<span class="p">[</span>g<span class="p">]</span><span class="o">+</span><span class="m">1</span>L<span class="p">)</span>:<span class="p">(</span>group.op.idx<span class="p">[</span>g<span class="o">+</span><span class="m">1</span><span class="p">]</span><span class="o">-</span><span class="m">1</span>L<span class="p">),]</span>    </div><div class='line' id='LC127'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;LIST.group <span class="o">&lt;-</span> getLIST<span class="p">(</span>FLAT.group<span class="p">,</span> meanstructure <span class="o">=</span> meanstructure<span class="p">,</span> </div><div class='line' id='LC128'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;int.ov.free <span class="o">=</span> int.ov.free<span class="p">,</span> int.lv.free <span class="o">=</span> int.lv.free<span class="p">,</span></div><div class='line' id='LC129'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;orthogonal <span class="o">=</span> orthogonal<span class="p">,</span> std.lv <span class="o">=</span> std.lv<span class="p">,</span> fixed.x <span class="o">=</span> fixed.x<span class="p">,</span></div><div class='line' id='LC130'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;auto.fix.first <span class="o">=</span> auto.fix.first<span class="p">,</span> </div><div class='line' id='LC131'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;auto.fix.single <span class="o">=</span> auto.fix.single<span class="p">,</span></div><div class='line' id='LC132'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;auto.var <span class="o">=</span> auto.var<span class="p">,</span> auto.cov.lv.x <span class="o">=</span> auto.cov.lv.x<span class="p">,</span></div><div class='line' id='LC133'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;auto.cov.y <span class="o">=</span> auto.cov.y<span class="p">,</span> group.equal <span class="o">=</span> <span class="kc">NULL</span><span class="p">,</span> ngroups <span class="o">=</span> <span class="m">1</span>L<span class="p">)</span></div><div class='line' id='LC134'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;LIST.group <span class="o">&lt;-</span> as.data.frame<span class="p">(</span>LIST.group<span class="p">,</span> stringsAsFactors<span class="o">=</span><span class="kc">FALSE</span><span class="p">)</span></div><div class='line' id='LC135'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>g <span class="o">==</span> <span class="m">1</span>L<span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC136'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;LIST <span class="o">&lt;-</span> LIST.group</div><div class='line' id='LC137'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span> <span class="kr">else</span> <span class="p">{</span></div><div class='line' id='LC138'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;LIST.group<span class="p">$</span>group <span class="o">&lt;-</span> rep<span class="p">(</span>g<span class="p">,</span> length<span class="p">(</span>LIST.group<span class="p">$</span>lhs<span class="p">))</span></div><div class='line' id='LC139'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;LIST.group<span class="p">$</span>id <span class="o">&lt;-</span> LIST.group<span class="p">$</span>id <span class="o">+</span> max<span class="p">(</span>LIST<span class="p">$</span>id<span class="p">)</span></div><div class='line' id='LC140'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;LIST <span class="o">&lt;-</span> rbind<span class="p">(</span>LIST<span class="p">,</span> LIST.group<span class="p">)</span></div><div class='line' id='LC141'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC142'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC143'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;LIST <span class="o">&lt;-</span> as.list<span class="p">(</span>LIST<span class="p">)</span></div><div class='line' id='LC144'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span> <span class="kr">else</span> <span class="p">{</span></div><div class='line' id='LC145'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;LIST <span class="o">&lt;-</span> getLIST<span class="p">(</span>FLAT<span class="p">,</span> meanstructure <span class="o">=</span> meanstructure<span class="p">,</span> </div><div class='line' id='LC146'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;int.ov.free <span class="o">=</span> int.ov.free<span class="p">,</span> int.lv.free <span class="o">=</span> int.lv.free<span class="p">,</span></div><div class='line' id='LC147'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;orthogonal <span class="o">=</span> orthogonal<span class="p">,</span> std.lv <span class="o">=</span> std.lv<span class="p">,</span> fixed.x <span class="o">=</span> fixed.x<span class="p">,</span></div><div class='line' id='LC148'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;auto.fix.first <span class="o">=</span> auto.fix.first<span class="p">,</span> auto.fix.single <span class="o">=</span> auto.fix.single<span class="p">,</span></div><div class='line' id='LC149'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;auto.var <span class="o">=</span> auto.var<span class="p">,</span> auto.cov.lv.x <span class="o">=</span> auto.cov.lv.x<span class="p">,</span></div><div class='line' id='LC150'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;auto.cov.y <span class="o">=</span> auto.cov.y<span class="p">,</span> group.equal <span class="o">=</span> group.equal<span class="p">,</span> </div><div class='line' id='LC151'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;ngroups <span class="o">=</span> ngroups<span class="p">)</span></div><div class='line' id='LC152'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span>        </div><div class='line' id='LC153'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>debug<span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC154'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;cat<span class="p">(</span><span class="s">&quot;[lavaan DEBUG]: parameter LIST:\n&quot;</span><span class="p">)</span></div><div class='line' id='LC155'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;print<span class="p">(</span> as.data.frame<span class="p">(</span>LIST<span class="p">,</span> stringsAsFactors<span class="o">=</span><span class="kc">FALSE</span><span class="p">)</span> <span class="p">)</span></div><div class='line' id='LC156'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC157'><br/></div><div class='line' id='LC158'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># apply user-specified modifiers</span></div><div class='line' id='LC159'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>length<span class="p">(</span>MOD<span class="p">))</span> <span class="p">{</span></div><div class='line' id='LC160'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">for</span><span class="p">(</span>el in <span class="m">1</span>:length<span class="p">(</span>MOD<span class="p">))</span> <span class="p">{</span></div><div class='line' id='LC161'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;idx <span class="o">&lt;-</span> which<span class="p">(</span>LIST<span class="p">$</span>mod.idx <span class="o">==</span> el<span class="p">)</span> <span class="c1"># for each group</span></div><div class='line' id='LC162'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;MOD.fixed <span class="o">&lt;-</span> MOD<span class="p">[[</span>el<span class="p">]]$</span>fixed  </div><div class='line' id='LC163'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;MOD.start <span class="o">&lt;-</span> MOD<span class="p">[[</span>el<span class="p">]]$</span>start</div><div class='line' id='LC164'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;MOD.label <span class="o">&lt;-</span> MOD<span class="p">[[</span>el<span class="p">]]$</span>label </div><div class='line' id='LC165'><br/></div><div class='line' id='LC166'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># check for single argument if multiple groups</span></div><div class='line' id='LC167'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>ngroups <span class="o">&gt;</span> <span class="m">1</span>L <span class="o">&amp;&amp;</span> length<span class="p">(</span>idx<span class="p">)</span> <span class="o">&gt;</span> <span class="m">1</span>L<span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC168'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># Ok, this is not very consistent:</span></div><div class='line' id='LC169'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># A) here we force same behavior across groups</span></div><div class='line' id='LC170'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>length<span class="p">(</span>MOD.fixed<span class="p">)</span> <span class="o">==</span> <span class="m">1</span>L<span class="p">)</span> MOD.fixed <span class="o">&lt;-</span> rep<span class="p">(</span>MOD.fixed<span class="p">,</span> ngroups<span class="p">)</span></div><div class='line' id='LC171'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>length<span class="p">(</span>MOD.start<span class="p">)</span> <span class="o">==</span> <span class="m">1</span>L<span class="p">)</span> MOD.start <span class="o">&lt;-</span> rep<span class="p">(</span>MOD.start<span class="p">,</span> ngroups<span class="p">)</span></div><div class='line' id='LC172'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># B) here we do NOT! otherwise, it would imply an equality </span></div><div class='line' id='LC173'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">#                    constraint...</span></div><div class='line' id='LC174'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">#    except if group.equal=&quot;loadings&quot;!</span></div><div class='line' id='LC175'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>length<span class="p">(</span>MOD.label<span class="p">)</span> <span class="o">==</span> <span class="m">1</span>L<span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC176'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span><span class="s">&quot;loadings&quot;</span> <span class="o">%in%</span> group.equal<span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC177'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;MOD.label <span class="o">&lt;-</span> rep<span class="p">(</span>MOD.label<span class="p">,</span> ngroups<span class="p">)</span></div><div class='line' id='LC178'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span> <span class="kr">else</span> <span class="p">{</span></div><div class='line' id='LC179'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;MOD.label <span class="o">&lt;-</span> c<span class="p">(</span>MOD.label<span class="p">,</span> rep<span class="p">(</span><span class="s">&quot;&quot;</span><span class="p">,</span> <span class="p">(</span>ngroups<span class="o">-</span><span class="m">1</span>L<span class="p">))</span> <span class="p">)</span></div><div class='line' id='LC180'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC181'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC182'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC183'><br/></div><div class='line' id='LC184'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># check for wrong number of arguments if multiple groups</span></div><div class='line' id='LC185'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;nidx <span class="o">&lt;-</span> length<span class="p">(</span>idx<span class="p">)</span></div><div class='line' id='LC186'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span> <span class="p">(</span>!is.null<span class="p">(</span>MOD.fixed<span class="p">)</span> <span class="o">&amp;&amp;</span> nidx <span class="o">!=</span> length<span class="p">(</span>MOD.fixed<span class="p">))</span> <span class="o">||</span></div><div class='line' id='LC187'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">(</span>!is.null<span class="p">(</span>MOD.start<span class="p">)</span> <span class="o">&amp;&amp;</span> nidx <span class="o">!=</span> length<span class="p">(</span>MOD.start<span class="p">))</span> <span class="o">||</span></div><div class='line' id='LC188'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">(</span>!is.null<span class="p">(</span>MOD.label<span class="p">)</span> <span class="o">&amp;&amp;</span> nidx <span class="o">!=</span> length<span class="p">(</span>MOD.label<span class="p">))</span> <span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC189'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;el.idx <span class="o">&lt;-</span> which<span class="p">(</span>LIST<span class="p">$</span>mod.idx <span class="o">==</span> el<span class="p">)[</span><span class="m">1</span>L<span class="p">]</span></div><div class='line' id='LC190'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;stop<span class="p">(</span><span class="s">&quot;lavaan ERROR: wrong number of arguments in modifier (&quot;</span><span class="p">,</span></div><div class='line' id='LC191'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;paste<span class="p">(</span>MOD.label<span class="p">,</span> collapse<span class="o">=</span><span class="s">&quot;,&quot;</span><span class="p">),</span> <span class="s">&quot;) of element &quot;</span><span class="p">,</span> </div><div class='line' id='LC192'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;LIST<span class="p">$</span>lhs<span class="p">[</span>el.idx<span class="p">],</span> LIST<span class="p">$</span>op<span class="p">[</span>el.idx<span class="p">],</span> LIST<span class="p">$</span>rhs<span class="p">[</span>el.idx<span class="p">])</span></div><div class='line' id='LC193'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC194'><br/></div><div class='line' id='LC195'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># apply modifiers</span></div><div class='line' id='LC196'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>!is.null<span class="p">(</span>MOD.fixed<span class="p">))</span> <span class="p">{</span></div><div class='line' id='LC197'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;na.idx <span class="o">&lt;-</span> which<span class="p">(</span>is.na<span class="p">(</span>MOD.fixed<span class="p">))</span></div><div class='line' id='LC198'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;not.na.idx <span class="o">&lt;-</span> which<span class="p">(</span>!is.na<span class="p">(</span>MOD.fixed<span class="p">))</span></div><div class='line' id='LC199'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;LIST<span class="p">$</span>ustart<span class="p">[</span>idx<span class="p">][</span>not.na.idx<span class="p">]</span> <span class="o">&lt;-</span> MOD.fixed<span class="p">[</span>not.na.idx<span class="p">]</span></div><div class='line' id='LC200'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;LIST<span class="p">$</span>free<span class="p">[</span>  idx<span class="p">][</span>not.na.idx<span class="p">]</span> <span class="o">&lt;-</span> <span class="m">0</span>L</div><div class='line' id='LC201'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;LIST<span class="p">$</span>free<span class="p">[</span>  idx<span class="p">][</span>na.idx<span class="p">]</span>     <span class="o">&lt;-</span> <span class="m">1</span>L <span class="c1"># eg factor loading</span></div><div class='line' id='LC202'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC203'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>!is.null<span class="p">(</span>MOD.start<span class="p">))</span> <span class="p">{</span></div><div class='line' id='LC204'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;LIST<span class="p">$</span>ustart<span class="p">[</span>idx<span class="p">]</span> <span class="o">&lt;-</span> MOD.start</div><div class='line' id='LC205'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC206'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>!is.null<span class="p">(</span>MOD.label<span class="p">))</span> <span class="p">{</span></div><div class='line' id='LC207'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;LIST<span class="p">$</span>label<span class="p">[</span>idx<span class="p">]</span> <span class="o">&lt;-</span> MOD.label</div><div class='line' id='LC208'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC209'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC210'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC211'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># remove mod.idx column</span></div><div class='line' id='LC212'>&nbsp;&nbsp;&nbsp;&nbsp;LIST<span class="p">$</span>mod.idx <span class="o">&lt;-</span> <span class="kc">NULL</span></div><div class='line' id='LC213'><br/></div><div class='line' id='LC214'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># check if CON contains *simple* equality constraints (eg b1 == b2)</span></div><div class='line' id='LC215'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># FIXME!!!</span></div><div class='line' id='LC216'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># b1 == b3</span></div><div class='line' id='LC217'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># b2 == b3 does not work</span></div><div class='line' id='LC218'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">#if(length(CON) &gt; 0L) {</span></div><div class='line' id='LC219'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">#    el.idx &lt;- integer(0L)</span></div><div class='line' id='LC220'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">#    for(el in 1:length(CON)) {</span></div><div class='line' id='LC221'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">#        if(CON[[el]]$op == &quot;==&quot;) {</span></div><div class='line' id='LC222'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">#            lhs.idx &lt;- which(LIST$label == CON[[el]]$lhs)</span></div><div class='line' id='LC223'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">#            rhs.idx &lt;- which(LIST$label == CON[[el]]$rhs)</span></div><div class='line' id='LC224'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">#            if(length(lhs.idx) &amp;&amp; length(rhs.idx)) {</span></div><div class='line' id='LC225'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">#                # flag this constraint (to be removed)</span></div><div class='line' id='LC226'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">#                el.idx &lt;- c(el.idx, el)</span></div><div class='line' id='LC227'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">#                # fill in equal and fix</span></div><div class='line' id='LC228'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">#                if(LIST$label[lhs.idx] == &quot;&quot;) {</span></div><div class='line' id='LC229'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">#                    LIST$label[rhs.idx] &lt;- LIST$label[lhs.idx]</span></div><div class='line' id='LC230'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">#                } else {</span></div><div class='line' id='LC231'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">#                    LIST$label[rhs.idx] &lt;- LIST$label[lhs.idx]</span></div><div class='line' id='LC232'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">#                }</span></div><div class='line' id='LC233'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">#                LIST$free[ rhs.idx] &lt;- 0L</span></div><div class='line' id='LC234'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">#            }</span></div><div class='line' id='LC235'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">#        }</span></div><div class='line' id='LC236'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">#    }</span></div><div class='line' id='LC237'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">#    if(length(el.idx) &gt; 0L) CON &lt;- CON[-el.idx]</span></div><div class='line' id='LC238'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">#}</span></div><div class='line' id='LC239'><br/></div><div class='line' id='LC240'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># get &#39;virtual&#39; parameter labels</span></div><div class='line' id='LC241'>&nbsp;&nbsp;&nbsp;&nbsp;LABEL <span class="o">&lt;-</span> getParameterLabels<span class="p">(</span>partable<span class="o">=</span>LIST<span class="p">,</span> group.equal<span class="o">=</span>group.equal<span class="p">,</span></div><div class='line' id='LC242'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;group.partial<span class="o">=</span>group.partial<span class="p">)</span></div><div class='line' id='LC243'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">#cat(&quot;DEBUG: label after getParameterLabels:\n&quot;); print(LABEL); cat(&quot;\n&quot;)</span></div><div class='line' id='LC244'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">#cat(&quot;DEBUG: eq.id after group.equal:\n&quot;); print(LIST$eq.id); cat(&quot;\n&quot;)</span></div><div class='line' id='LC245'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">#cat(&quot;DEBUG: LIST$label:\n&quot;); print(LIST$label); cat(&quot;\n&quot;)</span></div><div class='line' id='LC246'><br/></div><div class='line' id='LC247'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># handle user-specified equality constraints</span></div><div class='line' id='LC248'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># insert &#39;target/reference&#39; id&#39;s in eq.id/label columns</span></div><div class='line' id='LC249'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># so that the eq.id column reflects sets of equal parameters</span></div><div class='line' id='LC250'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># the &#39;reference&#39; parameter has free &gt; 0; the &#39;equal to&#39; pars free == 0</span></div><div class='line' id='LC251'>&nbsp;&nbsp;&nbsp;&nbsp;idx.eq.label <span class="o">&lt;-</span> which<span class="p">(</span>duplicated<span class="p">(</span>LABEL<span class="p">))</span></div><div class='line' id='LC252'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>length<span class="p">(</span>idx.eq.label<span class="p">)</span> <span class="o">&gt;</span> <span class="m">0</span>L<span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC253'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">for</span><span class="p">(</span>idx in idx.eq.label<span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC254'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;eq.label <span class="o">&lt;-</span> LABEL<span class="p">[</span>idx<span class="p">]</span></div><div class='line' id='LC255'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;ref.idx <span class="o">&lt;-</span> which<span class="p">(</span>LABEL <span class="o">==</span> eq.label<span class="p">)[</span><span class="m">1</span>L<span class="p">]</span> <span class="c1"># the first one only</span></div><div class='line' id='LC256'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># set eq.id equal</span></div><div class='line' id='LC257'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;LIST<span class="p">$</span>eq.id<span class="p">[</span>ref.idx<span class="p">]</span> <span class="o">&lt;-</span> LIST<span class="p">$</span>eq.id<span class="p">[</span>idx<span class="p">]</span> <span class="o">&lt;-</span> ref.idx</div><div class='line' id='LC258'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># fix target</span></div><div class='line' id='LC259'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;LIST<span class="p">$</span>free<span class="p">[</span>idx<span class="p">]</span> <span class="o">&lt;-</span> <span class="m">0</span>L</div><div class='line' id='LC260'><br/></div><div class='line' id='LC261'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># special case: check if there are any more instances </span></div><div class='line' id='LC262'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># of idx in  LIST$eq.id (perhaps due to group.equal)</span></div><div class='line' id='LC263'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;idx.all <span class="o">&lt;-</span> which<span class="p">(</span>LIST<span class="p">$</span>eq.id <span class="o">==</span> idx<span class="p">)</span></div><div class='line' id='LC264'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>length<span class="p">(</span>idx.all<span class="p">)</span> <span class="o">&gt;</span> <span class="m">0</span>L<span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC265'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;ref.idx.all <span class="o">&lt;-</span> which<span class="p">(</span>LIST<span class="p">$</span>eq.id <span class="o">==</span> ref.idx<span class="p">)</span></div><div class='line' id='LC266'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;LIST<span class="p">$</span>label<span class="p">[</span>ref.idx.all<span class="p">]</span> <span class="o">&lt;-</span> eq.label</div><div class='line' id='LC267'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;LIST<span class="p">$</span>eq.id<span class="p">[</span>idx.all<span class="p">]</span> <span class="o">&lt;-</span> ref.idx</div><div class='line' id='LC268'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;LIST<span class="p">$</span>free<span class="p">[</span>idx.all<span class="p">]</span> <span class="o">&lt;-</span> <span class="m">0</span>L  <span class="c1"># fix!</span></div><div class='line' id='LC269'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;LIST<span class="p">$</span>label<span class="p">[</span>idx.all<span class="p">]</span> <span class="o">&lt;-</span> eq.label</div><div class='line' id='LC270'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC271'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># special case: ref.idx is a fixed parameter</span></div><div class='line' id='LC272'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>LIST<span class="p">$</span>free<span class="p">[</span>ref.idx<span class="p">]</span> <span class="o">==</span> <span class="m">0</span>L<span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC273'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;ref.idx.all <span class="o">&lt;-</span> which<span class="p">(</span>LIST<span class="p">$</span>eq.id <span class="o">==</span> ref.idx<span class="p">)</span></div><div class='line' id='LC274'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;LIST<span class="p">$</span>ustart<span class="p">[</span>ref.idx.all<span class="p">]</span> <span class="o">&lt;-</span> LIST<span class="p">$</span>ustart<span class="p">[</span>ref.idx<span class="p">]</span></div><div class='line' id='LC275'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;LIST<span class="p">$</span>eq.id<span class="p">[</span>ref.idx.all<span class="p">]</span> <span class="o">&lt;-</span> <span class="m">0</span>L</div><div class='line' id='LC276'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC277'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC278'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC279'>&nbsp;</div><div class='line' id='LC280'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">#cat(&quot;DEBUG: eq.id after eq.id:\n&quot;); print(LIST$eq.id); cat(&quot;\n&quot;)    </span></div><div class='line' id='LC281'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">#cat(&quot;DEBUG: LIST$label:\n&quot;); print(LIST$label); cat(&quot;\n&quot;)</span></div><div class='line' id='LC282'><br/></div><div class='line' id='LC283'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># count free parameters</span></div><div class='line' id='LC284'>&nbsp;&nbsp;&nbsp;&nbsp;idx.free <span class="o">&lt;-</span> which<span class="p">(</span>LIST<span class="p">$</span>free <span class="o">&gt;</span> <span class="m">0</span><span class="p">)</span></div><div class='line' id='LC285'>&nbsp;&nbsp;&nbsp;&nbsp;LIST<span class="p">$</span>free<span class="p">[</span>idx.free<span class="p">]</span> <span class="o">&lt;-</span> <span class="m">1</span>:length<span class="p">(</span>idx.free<span class="p">)</span></div><div class='line' id='LC286'><br/></div><div class='line' id='LC287'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># 2. add free counter to this element</span></div><div class='line' id='LC288'>&nbsp;&nbsp;&nbsp;&nbsp;idx.equal <span class="o">&lt;-</span> which<span class="p">(</span>LIST<span class="p">$</span>eq.id <span class="o">&gt;</span> <span class="m">0</span><span class="p">)</span></div><div class='line' id='LC289'>&nbsp;&nbsp;&nbsp;&nbsp;LIST<span class="p">$</span>free<span class="p">[</span>idx.equal<span class="p">]</span> <span class="o">&lt;-</span> LIST<span class="p">$</span>free<span class="p">[</span> LIST<span class="p">$</span>eq.id<span class="p">[</span>idx.equal<span class="p">]</span> <span class="p">]</span></div><div class='line' id='LC290'><br/></div><div class='line' id='LC291'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># 3. which parameters would be free without equality constraints?</span></div><div class='line' id='LC292'>&nbsp;&nbsp;&nbsp;&nbsp;idx.unco <span class="o">&lt;-</span> which<span class="p">(</span>LIST<span class="p">$</span>free <span class="o">&gt;</span> <span class="m">0</span><span class="p">)</span></div><div class='line' id='LC293'>&nbsp;&nbsp;&nbsp;&nbsp;LIST<span class="p">$</span>unco<span class="p">[</span>idx.unco<span class="p">]</span> <span class="o">&lt;-</span> <span class="m">1</span>:length<span class="p">(</span>idx.unco<span class="p">)</span></div><div class='line' id='LC294'><br/></div><div class='line' id='LC295'><br/></div><div class='line' id='LC296'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># handle constraints (if any) (NOT per group, but overall - 0.4-11)</span></div><div class='line' id='LC297'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>length<span class="p">(</span>CON<span class="p">)</span> <span class="o">&gt;</span> <span class="m">0</span>L<span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC298'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">#cat(&quot;DEBUG:\n&quot;); print(CON)</span></div><div class='line' id='LC299'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;lhs <span class="o">=</span> unlist<span class="p">(</span>lapply<span class="p">(</span>CON<span class="p">,</span> <span class="s">&quot;[[&quot;</span><span class="p">,</span> <span class="s">&quot;lhs&quot;</span><span class="p">))</span></div><div class='line' id='LC300'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;op <span class="o">=</span> unlist<span class="p">(</span>lapply<span class="p">(</span>CON<span class="p">,</span> <span class="s">&quot;[[&quot;</span><span class="p">,</span>  <span class="s">&quot;op&quot;</span><span class="p">))</span></div><div class='line' id='LC301'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;rhs <span class="o">=</span> unlist<span class="p">(</span>lapply<span class="p">(</span>CON<span class="p">,</span> <span class="s">&quot;[[&quot;</span><span class="p">,</span> <span class="s">&quot;rhs&quot;</span><span class="p">))</span></div><div class='line' id='LC302'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;LIST<span class="p">$</span>id         <span class="o">&lt;-</span> c<span class="p">(</span>LIST<span class="p">$</span>id<span class="p">,</span>         max<span class="p">(</span>LIST<span class="p">$</span>id<span class="p">)</span> <span class="o">+</span> <span class="m">1</span>:length<span class="p">(</span>lhs<span class="p">)</span> <span class="p">)</span></div><div class='line' id='LC303'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;LIST<span class="p">$</span>lhs        <span class="o">&lt;-</span> c<span class="p">(</span>LIST<span class="p">$</span>lhs<span class="p">,</span>        lhs<span class="p">)</span></div><div class='line' id='LC304'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;LIST<span class="p">$</span>op         <span class="o">&lt;-</span> c<span class="p">(</span>LIST<span class="p">$</span>op<span class="p">,</span>         op<span class="p">)</span></div><div class='line' id='LC305'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;LIST<span class="p">$</span>rhs        <span class="o">&lt;-</span> c<span class="p">(</span>LIST<span class="p">$</span>rhs<span class="p">,</span>        rhs<span class="p">)</span></div><div class='line' id='LC306'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;LIST<span class="p">$</span>user       <span class="o">&lt;-</span> c<span class="p">(</span>LIST<span class="p">$</span>user<span class="p">,</span>       rep<span class="p">(</span><span class="m">1</span>L<span class="p">,</span> length<span class="p">(</span>lhs<span class="p">))</span> <span class="p">)</span></div><div class='line' id='LC307'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;LIST<span class="p">$</span>group      <span class="o">&lt;-</span> c<span class="p">(</span>LIST<span class="p">$</span>group<span class="p">,</span>      rep<span class="p">(</span><span class="m">0</span>L<span class="p">,</span> each<span class="o">=</span>length<span class="p">(</span>CON<span class="p">)))</span></div><div class='line' id='LC308'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;LIST<span class="p">$</span>free       <span class="o">&lt;-</span> c<span class="p">(</span>LIST<span class="p">$</span>free<span class="p">,</span>       rep<span class="p">(</span><span class="m">0</span>L<span class="p">,</span> length<span class="p">(</span>lhs<span class="p">))</span> <span class="p">)</span></div><div class='line' id='LC309'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;LIST<span class="p">$</span>ustart     <span class="o">&lt;-</span> c<span class="p">(</span>LIST<span class="p">$</span>ustart<span class="p">,</span>     rep<span class="p">(</span>as.numeric<span class="p">(</span><span class="kc">NA</span><span class="p">),</span> length<span class="p">(</span>lhs<span class="p">)))</span></div><div class='line' id='LC310'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;LIST<span class="p">$</span>exo        <span class="o">&lt;-</span> c<span class="p">(</span>LIST<span class="p">$</span>exo<span class="p">,</span>        rep<span class="p">(</span><span class="m">0</span>L<span class="p">,</span> length<span class="p">(</span>lhs<span class="p">))</span> <span class="p">)</span></div><div class='line' id='LC311'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;LIST<span class="p">$</span>label      <span class="o">&lt;-</span> c<span class="p">(</span>LIST<span class="p">$</span>label<span class="p">,</span>      rep<span class="p">(</span><span class="s">&quot;&quot;</span><span class="p">,</span>  length<span class="p">(</span>lhs<span class="p">))</span> <span class="p">)</span></div><div class='line' id='LC312'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;LIST<span class="p">$</span>eq.id      <span class="o">&lt;-</span> c<span class="p">(</span>LIST<span class="p">$</span>eq.id<span class="p">,</span>      rep<span class="p">(</span><span class="m">0</span>L<span class="p">,</span>  length<span class="p">(</span>lhs<span class="p">))</span> <span class="p">)</span></div><div class='line' id='LC313'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;LIST<span class="p">$</span>unco       <span class="o">&lt;-</span> c<span class="p">(</span>LIST<span class="p">$</span>unco<span class="p">,</span>       rep<span class="p">(</span><span class="m">0</span>L<span class="p">,</span>  length<span class="p">(</span>lhs<span class="p">))</span> <span class="p">)</span></div><div class='line' id='LC314'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC315'><br/></div><div class='line' id='LC316'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># put lhs of := elements in label column</span></div><div class='line' id='LC317'>&nbsp;&nbsp;&nbsp;&nbsp;def.idx <span class="o">&lt;-</span> which<span class="p">(</span>LIST<span class="p">$</span>op <span class="o">==</span> <span class="s">&quot;:=&quot;</span><span class="p">)</span></div><div class='line' id='LC318'>&nbsp;&nbsp;&nbsp;&nbsp;LIST<span class="p">$</span>label<span class="p">[</span>def.idx<span class="p">]</span> <span class="o">&lt;-</span> LIST<span class="p">$</span>lhs<span class="p">[</span>def.idx<span class="p">]</span></div><div class='line' id='LC319'><br/></div><div class='line' id='LC320'><br/></div><div class='line' id='LC321'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>debug<span class="p">)</span> <span class="p">{</span> </div><div class='line' id='LC322'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;cat<span class="p">(</span><span class="s">&quot;[lavaan DEBUG] lavParTable\n&quot;</span><span class="p">)</span></div><div class='line' id='LC323'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;print<span class="p">(</span> as.data.frame<span class="p">(</span>LIST<span class="p">)</span> <span class="p">)</span></div><div class='line' id='LC324'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC325'><br/></div><div class='line' id='LC326'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># data.frame?</span></div><div class='line' id='LC327'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>as.data.frame.<span class="p">)</span> LIST <span class="o">&lt;-</span> as.data.frame<span class="p">(</span>LIST<span class="p">,</span> stringsAsFactors <span class="o">=</span> <span class="kc">FALSE</span><span class="p">)</span></div><div class='line' id='LC328'><br/></div><div class='line' id='LC329'>&nbsp;&nbsp;&nbsp;&nbsp;LIST</div><div class='line' id='LC330'><span class="p">}</span></div><div class='line' id='LC331'><br/></div><div class='line' id='LC332'>parseModelString <span class="o">&lt;-</span> <span class="kr">function</span><span class="p">(</span>model.syntax <span class="o">=</span> <span class="s">&#39;&#39;</span><span class="p">,</span> as.data.frame. <span class="o">=</span> <span class="kc">FALSE</span><span class="p">,</span></div><div class='line' id='LC333'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;warn <span class="o">=</span> <span class="kc">TRUE</span><span class="p">,</span> debug <span class="o">=</span> <span class="kc">FALSE</span><span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC334'><br/></div><div class='line' id='LC335'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># check for empty syntax</span></div><div class='line' id='LC336'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>length<span class="p">(</span>model.syntax<span class="p">)</span> <span class="o">==</span> <span class="m">0</span><span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC337'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;stop<span class="p">(</span><span class="s">&quot;lavaan ERROR: empty model syntax&quot;</span><span class="p">)</span></div><div class='line' id='LC338'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC339'><br/></div><div class='line' id='LC340'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># break up in lines </span></div><div class='line' id='LC341'>&nbsp;&nbsp;&nbsp;&nbsp;model <span class="o">&lt;-</span> unlist<span class="p">(</span> strsplit<span class="p">(</span>model.syntax<span class="p">,</span> <span class="s">&quot;\n&quot;</span><span class="p">)</span> <span class="p">)</span></div><div class='line' id='LC342'><br/></div><div class='line' id='LC343'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># remove comments starting with &#39;#&#39; or &#39;!&#39;</span></div><div class='line' id='LC344'>&nbsp;&nbsp;&nbsp;&nbsp;model <span class="o">&lt;-</span> gsub<span class="p">(</span><span class="s">&quot;#.*&quot;</span><span class="p">,</span><span class="s">&quot;&quot;</span><span class="p">,</span> model<span class="p">);</span> model <span class="o">&lt;-</span> gsub<span class="p">(</span><span class="s">&quot;!.*&quot;</span><span class="p">,</span><span class="s">&quot;&quot;</span><span class="p">,</span> model<span class="p">)</span></div><div class='line' id='LC345'><br/></div><div class='line' id='LC346'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># replace semicolons by newlines and split in lines again</span></div><div class='line' id='LC347'>&nbsp;&nbsp;&nbsp;&nbsp;model <span class="o">&lt;-</span> gsub<span class="p">(</span><span class="s">&quot;;&quot;</span><span class="p">,</span><span class="s">&quot;\n&quot;</span><span class="p">,</span> model<span class="p">);</span> model <span class="o">&lt;-</span> unlist<span class="p">(</span> strsplit<span class="p">(</span>model<span class="p">,</span> <span class="s">&quot;\n&quot;</span><span class="p">)</span> <span class="p">)</span></div><div class='line' id='LC348'><br/></div><div class='line' id='LC349'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># strip all white space</span></div><div class='line' id='LC350'>&nbsp;&nbsp;&nbsp;&nbsp;model <span class="o">&lt;-</span> gsub<span class="p">(</span><span class="s">&quot;[[:space:]]+&quot;</span><span class="p">,</span> <span class="s">&quot;&quot;</span><span class="p">,</span> model<span class="p">)</span></div><div class='line' id='LC351'><br/></div><div class='line' id='LC352'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># keep non-empty lines only</span></div><div class='line' id='LC353'>&nbsp;&nbsp;&nbsp;&nbsp;idx <span class="o">&lt;-</span> which<span class="p">(</span>nzchar<span class="p">(</span>model<span class="p">))</span></div><div class='line' id='LC354'>&nbsp;&nbsp;&nbsp;&nbsp;model <span class="o">&lt;-</span> model<span class="p">[</span>idx<span class="p">]</span></div><div class='line' id='LC355'><br/></div><div class='line' id='LC356'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># check for multi-line formulas: they contain no &quot;~&quot; or &quot;=&quot; character</span></div><div class='line' id='LC357'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># but before we do that, we remove all modifiers</span></div><div class='line' id='LC358'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># to avoid confusion with for example equal(&quot;f1=~x1&quot;) statements</span></div><div class='line' id='LC359'>&nbsp;&nbsp;&nbsp;&nbsp;model.simple <span class="o">&lt;-</span> gsub<span class="p">(</span><span class="s">&quot;\\(.*\\)\\*&quot;</span><span class="p">,</span> <span class="s">&quot;MODIFIER*&quot;</span><span class="p">,</span> model<span class="p">)</span></div><div class='line' id='LC360'><br/></div><div class='line' id='LC361'>&nbsp;&nbsp;&nbsp;&nbsp;start.idx <span class="o">&lt;-</span> grep<span class="p">(</span><span class="s">&quot;[~=&lt;&gt;:]&quot;</span><span class="p">,</span> model.simple<span class="p">)</span></div><div class='line' id='LC362'>&nbsp;&nbsp;&nbsp;&nbsp;end.idx <span class="o">&lt;-</span> c<span class="p">(</span> start.idx<span class="p">[</span><span class="m">-1</span><span class="p">]</span><span class="o">-</span><span class="m">1</span><span class="p">,</span> length<span class="p">(</span>model<span class="p">)</span> <span class="p">)</span></div><div class='line' id='LC363'>&nbsp;&nbsp;&nbsp;&nbsp;model.orig    <span class="o">&lt;-</span> model</div><div class='line' id='LC364'>&nbsp;&nbsp;&nbsp;&nbsp;model <span class="o">&lt;-</span> character<span class="p">(</span> length<span class="p">(</span>start.idx<span class="p">)</span> <span class="p">)</span></div><div class='line' id='LC365'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">for</span><span class="p">(</span>i in <span class="m">1</span>:length<span class="p">(</span>start.idx<span class="p">))</span> <span class="p">{</span></div><div class='line' id='LC366'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;model<span class="p">[</span>i<span class="p">]</span> <span class="o">&lt;-</span> paste<span class="p">(</span>model.orig<span class="p">[</span>start.idx<span class="p">[</span>i<span class="p">]</span>:end.idx<span class="p">[</span>i<span class="p">]],</span> collapse<span class="o">=</span><span class="s">&quot;&quot;</span><span class="p">)</span></div><div class='line' id='LC367'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC368'><br/></div><div class='line' id='LC369'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># ok, in all remaining lines, we should have a &#39;~&#39; operator</span></div><div class='line' id='LC370'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># OR one of &#39;=&#39;, &#39;&lt;&#39; &#39;&gt;&#39; outside the &quot;&quot;</span></div><div class='line' id='LC371'>&nbsp;&nbsp;&nbsp;&nbsp;model.simple <span class="o">&lt;-</span> gsub<span class="p">(</span><span class="s">&quot;\\\&quot;</span><span class="m">.</span><span class="p">[</span><span class="o">^</span>\\\<span class="s">&quot;]*\\\&quot;&quot;, &quot;</span>LABEL<span class="s">&quot;, model)</span></div><div class='line' id='LC372'><span class="s">    idx.wrong &lt;- which(!grepl(&quot;</span><span class="p">[</span>~<span class="o">=&lt;&gt;</span>:<span class="p">]</span><span class="s">&quot;, model.simple))</span></div><div class='line' id='LC373'><span class="s">    if(length(idx.wrong) &gt; 0) {</span></div><div class='line' id='LC374'><span class="s">        cat(&quot;</span>lavaan: missing operator in formula<span class="p">(</span>s<span class="p">)</span>:\n<span class="s">&quot;)</span></div><div class='line' id='LC375'><span class="s">        print(model[idx.wrong])</span></div><div class='line' id='LC376'><span class="s">        stop(&quot;</span>lavaan ERROR: syntax error in lavaan model syntax<span class="s">&quot;)</span></div><div class='line' id='LC377'><span class="s">    }</span></div><div class='line' id='LC378'><span class="s">  </span></div><div class='line' id='LC379'><br/></div><div class='line' id='LC380'><span class="s">    # main operation: flatten formulas into single bivariate pieces</span></div><div class='line' id='LC381'><span class="s">    # with a left-hand-side (lhs), an operator (eg &quot;</span><span class="o">=</span>~<span class="s">&quot;), and a </span></div><div class='line' id='LC382'><span class="s">    # right-hand-side (rhs)</span></div><div class='line' id='LC383'><span class="s">    # both lhs and rhs can have a modifier </span></div><div class='line' id='LC384'><span class="s">    # (but we ignore the lhs modifier for now)</span></div><div class='line' id='LC385'><span class="s">    FLAT.lhs         &lt;- character(0)</span></div><div class='line' id='LC386'><span class="s">    #FLAT.lhs.mod    &lt;- character(0)</span></div><div class='line' id='LC387'><span class="s">    FLAT.op          &lt;- character(0)</span></div><div class='line' id='LC388'><span class="s">    FLAT.rhs         &lt;- character(0)</span></div><div class='line' id='LC389'><span class="s">    FLAT.rhs.mod.idx &lt;- integer(0)</span></div><div class='line' id='LC390'><span class="s">    FLAT.group       &lt;- integer(0)    # keep track of groups using &quot;</span>:<span class="s">&quot; operator</span></div><div class='line' id='LC391'><br/></div><div class='line' id='LC392'><span class="s">    FLAT.fixed       &lt;- character(0)  # only for display purposes! </span></div><div class='line' id='LC393'><span class="s">    FLAT.start       &lt;- character(0)  # only for display purposes!</span></div><div class='line' id='LC394'><span class="s">    FLAT.label       &lt;- character(0)  # only for display purposes!</span></div><div class='line' id='LC395'><span class="s">    FLAT.idx &lt;- 0L</span></div><div class='line' id='LC396'><span class="s">    MOD.idx  &lt;- 0L</span></div><div class='line' id='LC397'><span class="s">    CON.idx  &lt;- 0L</span></div><div class='line' id='LC398'><span class="s">    MOD &lt;- vector(&quot;</span>list<span class="s">&quot;, length=0L)</span></div><div class='line' id='LC399'><span class="s">    CON &lt;- vector(&quot;</span>list<span class="s">&quot;, length=0L)</span></div><div class='line' id='LC400'><span class="s">    GRP &lt;- 0L</span></div><div class='line' id='LC401'><span class="s">    for(i in 1:length(model)) {</span></div><div class='line' id='LC402'><span class="s">        x &lt;- model[i]</span></div><div class='line' id='LC403'><span class="s">        if(debug) {</span></div><div class='line' id='LC404'><span class="s">            cat(&quot;</span>formula to parse:\n<span class="s">&quot;); print(x); cat(&quot;</span>\n<span class="s">&quot;)</span></div><div class='line' id='LC405'><span class="s">        }</span></div><div class='line' id='LC406'><br/></div><div class='line' id='LC407'><span class="s">        # 1. which operator is used?</span></div><div class='line' id='LC408'><span class="s">        line.simple &lt;- gsub(&quot;</span>\\\<span class="s">&quot;.[^\\\&quot;</span><span class="p">]</span><span class="o">*</span>\\\<span class="s">&quot;&quot;</span><span class="p">,</span> <span class="s">&quot;LABEL&quot;</span><span class="p">,</span> x<span class="p">)</span></div><div class='line' id='LC409'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># &quot;=~&quot; operator?</span></div><div class='line' id='LC410'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>grepl<span class="p">(</span><span class="s">&quot;=~&quot;</span><span class="p">,</span> line.simple<span class="p">,</span> fixed<span class="o">=</span><span class="kc">TRUE</span><span class="p">))</span> <span class="p">{</span></div><div class='line' id='LC411'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;op <span class="o">&lt;-</span> <span class="s">&quot;=~&quot;</span></div><div class='line' id='LC412'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># &quot;&lt;~&quot; operator?</span></div><div class='line' id='LC413'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span> <span class="kr">else</span> <span class="kr">if</span><span class="p">(</span>grepl<span class="p">(</span><span class="s">&quot;&lt;~&quot;</span><span class="p">,</span> line.simple<span class="p">,</span> fixed<span class="o">=</span><span class="kc">TRUE</span><span class="p">))</span> <span class="p">{</span></div><div class='line' id='LC414'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;op <span class="o">&lt;-</span> <span class="s">&quot;&lt;~&quot;</span></div><div class='line' id='LC415'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># &quot;~~&quot; operator?</span></div><div class='line' id='LC416'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span> <span class="kr">else</span> <span class="kr">if</span><span class="p">(</span>grepl<span class="p">(</span><span class="s">&quot;~~&quot;</span><span class="p">,</span> line.simple<span class="p">,</span> fixed<span class="o">=</span><span class="kc">TRUE</span><span class="p">))</span> <span class="p">{</span></div><div class='line' id='LC417'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;op <span class="o">&lt;-</span> <span class="s">&quot;~~&quot;</span></div><div class='line' id='LC418'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># &quot;~&quot; operator?</span></div><div class='line' id='LC419'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span> <span class="kr">else</span> <span class="kr">if</span><span class="p">(</span>grepl<span class="p">(</span><span class="s">&quot;~&quot;</span><span class="p">,</span> line.simple<span class="p">,</span> fixed<span class="o">=</span><span class="kc">TRUE</span><span class="p">))</span> <span class="p">{</span></div><div class='line' id='LC420'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;op <span class="o">&lt;-</span> <span class="s">&quot;~&quot;</span>           </div><div class='line' id='LC421'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># &quot;==&quot; operator?</span></div><div class='line' id='LC422'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span> <span class="kr">else</span> <span class="kr">if</span><span class="p">(</span>grepl<span class="p">(</span><span class="s">&quot;==&quot;</span><span class="p">,</span> line.simple<span class="p">,</span> fixed<span class="o">=</span><span class="kc">TRUE</span><span class="p">))</span> <span class="p">{</span></div><div class='line' id='LC423'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;op <span class="o">&lt;-</span> <span class="s">&quot;==&quot;</span>  </div><div class='line' id='LC424'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># &quot;&lt;&quot; operator?</span></div><div class='line' id='LC425'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span> <span class="kr">else</span> <span class="kr">if</span><span class="p">(</span>grepl<span class="p">(</span><span class="s">&quot;&lt;&quot;</span><span class="p">,</span> line.simple<span class="p">,</span> fixed<span class="o">=</span><span class="kc">TRUE</span><span class="p">))</span> <span class="p">{</span></div><div class='line' id='LC426'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;op <span class="o">&lt;-</span> <span class="s">&quot;&lt;&quot;</span></div><div class='line' id='LC427'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># &quot;&gt;&quot; operator?</span></div><div class='line' id='LC428'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span> <span class="kr">else</span> <span class="kr">if</span><span class="p">(</span>grepl<span class="p">(</span><span class="s">&quot;&gt;&quot;</span><span class="p">,</span> line.simple<span class="p">,</span> fixed<span class="o">=</span><span class="kc">TRUE</span><span class="p">))</span> <span class="p">{</span></div><div class='line' id='LC429'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;op <span class="o">&lt;-</span> <span class="s">&quot;&gt;&quot;</span></div><div class='line' id='LC430'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># &quot;:=&quot; operator?</span></div><div class='line' id='LC431'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span> <span class="kr">else</span> <span class="kr">if</span><span class="p">(</span>grepl<span class="p">(</span><span class="s">&quot;:=&quot;</span><span class="p">,</span> line.simple<span class="p">,</span> fixed<span class="o">=</span><span class="kc">TRUE</span><span class="p">))</span> <span class="p">{</span></div><div class='line' id='LC432'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;op <span class="o">&lt;-</span> <span class="s">&quot;:=&quot;</span></div><div class='line' id='LC433'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># &quot;:&quot; operator?</span></div><div class='line' id='LC434'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span> <span class="kr">else</span> <span class="kr">if</span><span class="p">(</span>grepl<span class="p">(</span><span class="s">&quot;:&quot;</span><span class="p">,</span> line.simple<span class="p">,</span> fixed<span class="o">=</span><span class="kc">TRUE</span><span class="p">))</span> <span class="p">{</span></div><div class='line' id='LC435'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;op <span class="o">&lt;-</span> <span class="s">&quot;:&quot;</span></div><div class='line' id='LC436'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span> <span class="kr">else</span> <span class="p">{</span></div><div class='line' id='LC437'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;stop<span class="p">(</span><span class="s">&quot;unknown operator in &quot;</span><span class="p">,</span> model<span class="p">[</span>i<span class="p">])</span></div><div class='line' id='LC438'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC439'><br/></div><div class='line' id='LC440'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># 2. split by operator (only the *first* occurence!)</span></div><div class='line' id='LC441'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># check first if equal/label modifier has been used on the LEFT!</span></div><div class='line' id='LC442'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>substr<span class="p">(</span>x<span class="p">,</span><span class="m">1</span><span class="p">,</span><span class="m">5</span><span class="p">)</span> <span class="o">==</span> <span class="s">&quot;label&quot;</span><span class="p">)</span> stop<span class="p">(</span><span class="s">&quot;label modifier can not be used on the left-hand side of the operator&quot;</span><span class="p">)</span></div><div class='line' id='LC443'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;op.idx <span class="o">&lt;-</span> regexpr<span class="p">(</span>op<span class="p">,</span> x<span class="p">)</span></div><div class='line' id='LC444'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;lhs <span class="o">&lt;-</span> substr<span class="p">(</span>x<span class="p">,</span> <span class="m">1</span>L<span class="p">,</span> op.idx<span class="o">-</span><span class="m">1</span>L<span class="p">)</span></div><div class='line' id='LC445'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;rhs <span class="o">&lt;-</span> substr<span class="p">(</span>x<span class="p">,</span> op.idx<span class="o">+</span>attr<span class="p">(</span>op.idx<span class="p">,</span> <span class="s">&quot;match.length&quot;</span><span class="p">),</span> nchar<span class="p">(</span>x<span class="p">))</span></div><div class='line' id='LC446'><br/></div><div class='line' id='LC447'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># 2b. if operator is &quot;==&quot; or &quot;&lt;&quot; or &quot;&gt;&quot; or &quot;:=&quot;, put it in CON</span></div><div class='line' id='LC448'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>op <span class="o">==</span> <span class="s">&quot;==&quot;</span> <span class="o">||</span> op <span class="o">==</span> <span class="s">&quot;&lt;&quot;</span> <span class="o">||</span> op <span class="o">==</span> <span class="s">&quot;&gt;&quot;</span> <span class="o">||</span> op <span class="o">==</span> <span class="s">&quot;:=&quot;</span><span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC449'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># remove quotes, if any</span></div><div class='line' id='LC450'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;lhs <span class="o">&lt;-</span> gsub<span class="p">(</span><span class="s">&quot;\\\&quot;&quot;, &quot;&quot;, lhs)</span></div><div class='line' id='LC451'><span class="s">            rhs &lt;- gsub(&quot;</span>\\\<span class="s">&quot;&quot;</span><span class="p">,</span> <span class="s">&quot;&quot;</span><span class="p">,</span> rhs<span class="p">)</span></div><div class='line' id='LC452'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;CON.idx <span class="o">&lt;-</span> CON.idx <span class="o">+</span> <span class="m">1</span>L</div><div class='line' id='LC453'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;CON<span class="p">[[</span>CON.idx<span class="p">]]</span> <span class="o">&lt;-</span> list<span class="p">(</span>op<span class="o">=</span>op<span class="p">,</span> lhs<span class="o">=</span>lhs<span class="p">,</span> rhs<span class="o">=</span>rhs<span class="p">)</span></div><div class='line' id='LC454'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;next</div><div class='line' id='LC455'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC456'><br/></div><div class='line' id='LC457'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># 2c if operator is &quot;:&quot;, put it in GRP</span></div><div class='line' id='LC458'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>op <span class="o">==</span> <span class="s">&quot;:&quot;</span><span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC459'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;FLAT.idx <span class="o">&lt;-</span> FLAT.idx <span class="o">+</span> <span class="m">1</span>L</div><div class='line' id='LC460'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;GRP <span class="o">&lt;-</span> GRP <span class="o">+</span> <span class="m">1</span>L</div><div class='line' id='LC461'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;FLAT.lhs<span class="p">[</span>FLAT.idx<span class="p">]</span> <span class="o">&lt;-</span> lhs</div><div class='line' id='LC462'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;FLAT.op<span class="p">[</span> FLAT.idx<span class="p">]</span> <span class="o">&lt;-</span> op</div><div class='line' id='LC463'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;FLAT.rhs<span class="p">[</span>FLAT.idx<span class="p">]</span> <span class="o">&lt;-</span> <span class="s">&quot;&quot;</span></div><div class='line' id='LC464'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;FLAT.fixed<span class="p">[</span>FLAT.idx<span class="p">]</span> <span class="o">&lt;-</span> <span class="s">&quot;&quot;</span></div><div class='line' id='LC465'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;FLAT.start<span class="p">[</span>FLAT.idx<span class="p">]</span> <span class="o">&lt;-</span> <span class="s">&quot;&quot;</span></div><div class='line' id='LC466'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;FLAT.label<span class="p">[</span>FLAT.idx<span class="p">]</span> <span class="o">&lt;-</span> <span class="s">&quot;&quot;</span></div><div class='line' id='LC467'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;FLAT.rhs.mod.idx<span class="p">[</span>FLAT.idx<span class="p">]</span> <span class="o">&lt;-</span> <span class="m">0</span>L</div><div class='line' id='LC468'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;FLAT.group<span class="p">[</span>FLAT.idx<span class="p">]</span> <span class="o">&lt;-</span> GRP</div><div class='line' id='LC469'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;next</div><div class='line' id='LC470'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC471'><br/></div><div class='line' id='LC472'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># 3. parse left hand</span></div><div class='line' id='LC473'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">#    lhs modifiers will be ignored for now</span></div><div class='line' id='LC474'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;lhs.formula <span class="o">&lt;-</span> as.formula<span class="p">(</span>paste<span class="p">(</span><span class="s">&quot;~&quot;</span><span class="p">,</span>lhs<span class="p">))</span></div><div class='line' id='LC475'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;out <span class="o">&lt;-</span> parse.rhs<span class="p">(</span>rhs<span class="o">=</span>lhs.formula<span class="p">[[</span><span class="m">2</span>L<span class="p">]])</span></div><div class='line' id='LC476'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;lhs.names <span class="o">&lt;-</span> names<span class="p">(</span>out<span class="p">)</span></div><div class='line' id='LC477'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># check if we have modifiers</span></div><div class='line' id='LC478'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>sum<span class="p">(</span>sapply<span class="p">(</span>out<span class="p">,</span> length<span class="p">))</span> <span class="o">&gt;</span> <span class="m">0</span>L<span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC479'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;warning<span class="p">(</span><span class="s">&quot;lavaan WARNING: left-hand side of formula below contains modifier:\n&quot;</span><span class="p">,</span> x<span class="p">,</span><span class="s">&quot;\n&quot;</span><span class="p">)</span></div><div class='line' id='LC480'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC481'><br/></div><div class='line' id='LC482'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># 4. parse rhs (as rhs of a single-sided formula)</span></div><div class='line' id='LC483'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;rhs.formula <span class="o">&lt;-</span> as.formula<span class="p">(</span>paste<span class="p">(</span><span class="s">&quot;~&quot;</span><span class="p">,</span>rhs<span class="p">))</span></div><div class='line' id='LC484'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;out <span class="o">&lt;-</span> parse.rhs<span class="p">(</span>rhs<span class="o">=</span>rhs.formula<span class="p">[[</span><span class="m">2</span>L<span class="p">]])</span></div><div class='line' id='LC485'><br/></div><div class='line' id='LC486'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># for each lhs element</span></div><div class='line' id='LC487'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">for</span><span class="p">(</span>l in <span class="m">1</span>:length<span class="p">(</span>lhs.names<span class="p">))</span> <span class="p">{</span></div><div class='line' id='LC488'><br/></div><div class='line' id='LC489'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># for each rhs element</span></div><div class='line' id='LC490'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">for</span><span class="p">(</span>j in <span class="m">1</span>:length<span class="p">(</span>out<span class="p">))</span> <span class="p">{</span></div><div class='line' id='LC491'><br/></div><div class='line' id='LC492'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># catch intercepts</span></div><div class='line' id='LC493'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">((</span>op <span class="o">==</span> <span class="s">&quot;~&quot;</span> <span class="o">||</span> op <span class="o">==</span> <span class="s">&quot;~1&quot;</span><span class="p">)</span> <span class="o">&amp;&amp;</span> names<span class="p">(</span>out<span class="p">)[</span>j<span class="p">]</span> <span class="o">==</span> <span class="s">&quot;intercept&quot;</span><span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC494'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;op <span class="o">&lt;-</span> <span class="s">&quot;~1&quot;</span></div><div class='line' id='LC495'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;rhs.name <span class="o">&lt;-</span> <span class="s">&quot;&quot;</span></div><div class='line' id='LC496'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span> <span class="kr">else</span> <span class="p">{</span></div><div class='line' id='LC497'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;rhs.name <span class="o">&lt;-</span> names<span class="p">(</span>out<span class="p">)[</span>j<span class="p">]</span></div><div class='line' id='LC498'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC499'><br/></div><div class='line' id='LC500'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># check if we not already have this combination (in this group)</span></div><div class='line' id='LC501'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># 1. asymmetric (=~, ~, ~1)</span></div><div class='line' id='LC502'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;idx <span class="o">&lt;-</span> which<span class="p">(</span>FLAT.lhs <span class="o">==</span> lhs.names<span class="p">[</span>l<span class="p">]</span> <span class="o">&amp;</span></div><div class='line' id='LC503'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;FLAT.op  <span class="o">==</span> op <span class="o">&amp;</span></div><div class='line' id='LC504'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;FLAT.group <span class="o">==</span> GRP <span class="o">&amp;</span></div><div class='line' id='LC505'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;FLAT.rhs <span class="o">==</span> rhs.name<span class="p">)</span></div><div class='line' id='LC506'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>length<span class="p">(</span>idx<span class="p">)</span> <span class="o">&gt;</span> <span class="m">0</span>L<span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC507'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;stop<span class="p">(</span><span class="s">&quot;lavaan ERROR: duplicate model element in: &quot;</span><span class="p">,</span> model<span class="p">[</span>i<span class="p">])</span></div><div class='line' id='LC508'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC509'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># 2. symmetric (~~)</span></div><div class='line' id='LC510'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;idx <span class="o">&lt;-</span> which<span class="p">(</span>FLAT.lhs <span class="o">==</span> rhs.name <span class="o">&amp;</span></div><div class='line' id='LC511'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;FLAT.op  <span class="o">==</span> <span class="s">&quot;~~&quot;</span> <span class="o">&amp;</span></div><div class='line' id='LC512'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;FLAT.group <span class="o">==</span> GRP <span class="o">&amp;</span></div><div class='line' id='LC513'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;FLAT.rhs <span class="o">==</span> lhs.names<span class="p">[</span>l<span class="p">])</span></div><div class='line' id='LC514'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>length<span class="p">(</span>idx<span class="p">)</span> <span class="o">&gt;</span> <span class="m">0</span>L<span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC515'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;stop<span class="p">(</span><span class="s">&quot;lavaan ERROR: duplicate model element in: &quot;</span><span class="p">,</span> model<span class="p">[</span>i<span class="p">])</span></div><div class='line' id='LC516'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC517'><br/></div><div class='line' id='LC518'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;FLAT.idx <span class="o">&lt;-</span> FLAT.idx <span class="o">+</span> <span class="m">1</span>L</div><div class='line' id='LC519'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;FLAT.lhs<span class="p">[</span>FLAT.idx<span class="p">]</span> <span class="o">&lt;-</span> lhs.names<span class="p">[</span>l<span class="p">]</span></div><div class='line' id='LC520'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;FLAT.op<span class="p">[</span> FLAT.idx<span class="p">]</span> <span class="o">&lt;-</span> op</div><div class='line' id='LC521'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;FLAT.rhs<span class="p">[</span>FLAT.idx<span class="p">]</span> <span class="o">&lt;-</span> rhs.name</div><div class='line' id='LC522'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;FLAT.group<span class="p">[</span>FLAT.idx<span class="p">]</span> <span class="o">&lt;-</span> GRP</div><div class='line' id='LC523'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;FLAT.fixed<span class="p">[</span>FLAT.idx<span class="p">]</span> <span class="o">&lt;-</span> <span class="s">&quot;&quot;</span></div><div class='line' id='LC524'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;FLAT.start<span class="p">[</span>FLAT.idx<span class="p">]</span> <span class="o">&lt;-</span> <span class="s">&quot;&quot;</span></div><div class='line' id='LC525'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;FLAT.label<span class="p">[</span>FLAT.idx<span class="p">]</span> <span class="o">&lt;-</span> <span class="s">&quot;&quot;</span></div><div class='line' id='LC526'><br/></div><div class='line' id='LC527'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;mod <span class="o">&lt;-</span> list<span class="p">()</span></div><div class='line' id='LC528'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;rhs.mod <span class="o">&lt;-</span> <span class="m">0</span>L</div><div class='line' id='LC529'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>length<span class="p">(</span>out<span class="p">[[</span>j<span class="p">]]$</span>fixed<span class="p">)</span> <span class="o">&gt;</span> <span class="m">0</span>L<span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC530'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;mod<span class="p">$</span>fixed <span class="o">&lt;-</span> out<span class="p">[[</span>j<span class="p">]]$</span>fixed</div><div class='line' id='LC531'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;FLAT.fixed<span class="p">[</span>FLAT.idx<span class="p">]</span> <span class="o">&lt;-</span> paste<span class="p">(</span>mod<span class="p">$</span>fixed<span class="p">,</span> collapse<span class="o">=</span><span class="s">&quot;;&quot;</span><span class="p">)</span></div><div class='line' id='LC532'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;rhs.mod <span class="o">&lt;-</span> <span class="m">1</span>L</div><div class='line' id='LC533'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC534'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>length<span class="p">(</span>out<span class="p">[[</span>j<span class="p">]]$</span>start<span class="p">)</span> <span class="o">&gt;</span> <span class="m">0</span>L<span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC535'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;mod<span class="p">$</span>start <span class="o">&lt;-</span> out<span class="p">[[</span>j<span class="p">]]$</span>start</div><div class='line' id='LC536'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;FLAT.start<span class="p">[</span>FLAT.idx<span class="p">]</span> <span class="o">&lt;-</span> paste<span class="p">(</span>mod<span class="p">$</span>start<span class="p">,</span> collapse<span class="o">=</span><span class="s">&quot;;&quot;</span><span class="p">)</span></div><div class='line' id='LC537'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;rhs.mod <span class="o">&lt;-</span> <span class="m">1</span>L</div><div class='line' id='LC538'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC539'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>length<span class="p">(</span>out<span class="p">[[</span>j<span class="p">]]$</span>label<span class="p">)</span> <span class="o">&gt;</span> <span class="m">0</span>L<span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC540'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;mod<span class="p">$</span>label <span class="o">&lt;-</span> out<span class="p">[[</span>j<span class="p">]]$</span>label</div><div class='line' id='LC541'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;FLAT.label<span class="p">[</span>FLAT.idx<span class="p">]</span> <span class="o">&lt;-</span> paste<span class="p">(</span>mod<span class="p">$</span>label<span class="p">,</span> collapse<span class="o">=</span><span class="s">&quot;;&quot;</span><span class="p">)</span></div><div class='line' id='LC542'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;rhs.mod <span class="o">&lt;-</span> <span class="m">1</span>L</div><div class='line' id='LC543'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC544'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>op <span class="o">==</span> <span class="s">&quot;~1&quot;</span> <span class="o">&amp;&amp;</span> rhs <span class="o">==</span> <span class="s">&quot;0&quot;</span><span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC545'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;mod<span class="p">$</span>fixed <span class="o">&lt;-</span> <span class="m">0</span></div><div class='line' id='LC546'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;FLAT.fixed<span class="p">[</span>FLAT.idx<span class="p">]</span> <span class="o">&lt;-</span> paste<span class="p">(</span>mod<span class="p">$</span>fixed<span class="p">,</span> collapse<span class="o">=</span><span class="s">&quot;;&quot;</span><span class="p">)</span></div><div class='line' id='LC547'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;rhs.mod <span class="o">&lt;-</span> <span class="m">1</span>L</div><div class='line' id='LC548'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC549'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>op <span class="o">==</span> <span class="s">&quot;=~&quot;</span> <span class="o">&amp;&amp;</span> rhs <span class="o">==</span> <span class="s">&quot;0&quot;</span><span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC550'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;mod<span class="p">$</span>fixed <span class="o">&lt;-</span> <span class="m">0</span></div><div class='line' id='LC551'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;FLAT.rhs<span class="p">[</span>FLAT.idx<span class="p">]</span> <span class="o">&lt;-</span> FLAT.lhs<span class="p">[</span>FLAT.idx<span class="p">]</span></div><div class='line' id='LC552'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;FLAT.fixed<span class="p">[</span>FLAT.idx<span class="p">]</span> <span class="o">&lt;-</span> paste<span class="p">(</span>mod<span class="p">$</span>fixed<span class="p">,</span> collapse<span class="o">=</span><span class="s">&quot;;&quot;</span><span class="p">)</span></div><div class='line' id='LC553'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;rhs.mod <span class="o">&lt;-</span> <span class="m">1</span>L</div><div class='line' id='LC554'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC555'><br/></div><div class='line' id='LC556'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;FLAT.rhs.mod.idx<span class="p">[</span>FLAT.idx<span class="p">]</span> <span class="o">&lt;-</span> rhs.mod</div><div class='line' id='LC557'><br/></div><div class='line' id='LC558'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>rhs.mod <span class="o">&gt;</span> <span class="m">0</span>L<span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC559'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;MOD.idx <span class="o">&lt;-</span> MOD.idx <span class="o">+</span> <span class="m">1</span>L</div><div class='line' id='LC560'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;MOD<span class="p">[[</span>MOD.idx<span class="p">]]</span> <span class="o">&lt;-</span> mod</div><div class='line' id='LC561'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC562'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span> <span class="c1"># rhs elements</span></div><div class='line' id='LC563'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span> <span class="c1"># lhs elements</span></div><div class='line' id='LC564'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span> <span class="c1"># model elements</span></div><div class='line' id='LC565'><br/></div><div class='line' id='LC566'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># enumerate modifier indices</span></div><div class='line' id='LC567'>&nbsp;&nbsp;&nbsp;&nbsp;mod.idx <span class="o">&lt;-</span> which<span class="p">(</span>FLAT.rhs.mod.idx <span class="o">&gt;</span> <span class="m">0</span>L<span class="p">)</span></div><div class='line' id='LC568'>&nbsp;&nbsp;&nbsp;&nbsp;FLAT.rhs.mod.idx<span class="p">[</span> mod.idx <span class="p">]</span> <span class="o">&lt;-</span> <span class="m">1</span>:length<span class="p">(</span>mod.idx<span class="p">)</span></div><div class='line' id='LC569'><br/></div><div class='line' id='LC570'>&nbsp;&nbsp;&nbsp;&nbsp;FLAT <span class="o">&lt;-</span> list<span class="p">(</span>lhs<span class="o">=</span>FLAT.lhs<span class="p">,</span> op<span class="o">=</span>FLAT.op<span class="p">,</span> rhs<span class="o">=</span>FLAT.rhs<span class="p">,</span></div><div class='line' id='LC571'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;mod.idx<span class="o">=</span>FLAT.rhs.mod.idx<span class="p">,</span> group<span class="o">=</span>FLAT.group<span class="p">,</span></div><div class='line' id='LC572'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;fixed<span class="o">=</span>FLAT.fixed<span class="p">,</span> start<span class="o">=</span>FLAT.start<span class="p">,</span></div><div class='line' id='LC573'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;label<span class="o">=</span>FLAT.label<span class="p">)</span></div><div class='line' id='LC574'><br/></div><div class='line' id='LC575'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>as.data.frame.<span class="p">)</span> FLAT <span class="o">&lt;-</span> as.data.frame<span class="p">(</span>FLAT<span class="p">,</span> stringsAsFactors<span class="o">=</span><span class="kc">FALSE</span><span class="p">)</span></div><div class='line' id='LC576'><br/></div><div class='line' id='LC577'>&nbsp;&nbsp;&nbsp;&nbsp;attr<span class="p">(</span>FLAT<span class="p">,</span> <span class="s">&quot;modifiers&quot;</span><span class="p">)</span> <span class="o">&lt;-</span> MOD</div><div class='line' id='LC578'>&nbsp;&nbsp;&nbsp;&nbsp;attr<span class="p">(</span>FLAT<span class="p">,</span> <span class="s">&quot;constraints&quot;</span><span class="p">)</span> <span class="o">&lt;-</span> CON</div><div class='line' id='LC579'><br/></div><div class='line' id='LC580'>&nbsp;&nbsp;&nbsp;&nbsp;FLAT</div><div class='line' id='LC581'><span class="p">}</span></div><div class='line' id='LC582'><br/></div><div class='line' id='LC583'>parse.rhs <span class="o">&lt;-</span> <span class="kr">function</span><span class="p">(</span>rhs<span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC584'><br/></div><div class='line' id='LC585'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># new version YR 15 dec 2011!</span></div><div class='line' id='LC586'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># - no &#39;equal&#39; field anymore (only labels!)</span></div><div class='line' id='LC587'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># - every modifier is evaluated</span></div><div class='line' id='LC588'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># - unquoted labels are allowed (eg. x1 + x2 + c(v1,v2,v3)*x3)</span></div><div class='line' id='LC589'><br/></div><div class='line' id='LC590'>&nbsp;&nbsp;&nbsp;&nbsp;getModifier <span class="o">&lt;-</span> <span class="kr">function</span><span class="p">(</span>mod<span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC591'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>length<span class="p">(</span>mod<span class="p">)</span> <span class="o">==</span> <span class="m">1</span>L<span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC592'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># three possibilites: 1) numeric, 2) NA, or 3) quoted character</span></div><div class='line' id='LC593'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span> is.numeric<span class="p">(</span>mod<span class="p">)</span> <span class="p">)</span> </div><div class='line' id='LC594'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">return</span><span class="p">(</span> list<span class="p">(</span>fixed<span class="o">=</span>mod<span class="p">)</span> <span class="p">)</span></div><div class='line' id='LC595'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span> is.na<span class="p">(</span>mod<span class="p">)</span> <span class="p">)</span> </div><div class='line' id='LC596'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">return</span><span class="p">(</span> list<span class="p">(</span>fixed<span class="o">=</span>as.numeric<span class="p">(</span><span class="kc">NA</span><span class="p">))</span> <span class="p">)</span></div><div class='line' id='LC597'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span> is.character<span class="p">(</span>mod<span class="p">)</span> <span class="p">)</span></div><div class='line' id='LC598'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">return</span><span class="p">(</span> list<span class="p">(</span>label<span class="o">=</span>mod<span class="p">)</span> <span class="p">)</span></div><div class='line' id='LC599'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span> <span class="kr">else</span> <span class="kr">if</span><span class="p">(</span>mod<span class="p">[[</span><span class="m">1</span>L<span class="p">]]</span> <span class="o">==</span> <span class="s">&quot;start&quot;</span><span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC600'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;cof <span class="o">&lt;-</span> unlist<span class="p">(</span>lapply<span class="p">(</span>as.list<span class="p">(</span>mod<span class="p">)[</span><span class="m">-1</span><span class="p">],</span> </div><div class='line' id='LC601'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;eval<span class="p">,</span> envir<span class="o">=</span><span class="kc">NULL</span><span class="p">,</span> enclos<span class="o">=</span><span class="kc">NULL</span><span class="p">))</span></div><div class='line' id='LC602'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">return</span><span class="p">(</span> list<span class="p">(</span>start<span class="o">=</span>cof<span class="p">)</span> <span class="p">)</span></div><div class='line' id='LC603'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span> <span class="kr">else</span> <span class="kr">if</span><span class="p">(</span>mod<span class="p">[[</span><span class="m">1</span>L<span class="p">]]</span> <span class="o">==</span> <span class="s">&quot;equal&quot;</span><span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC604'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;label <span class="o">&lt;-</span> unlist<span class="p">(</span>lapply<span class="p">(</span>as.list<span class="p">(</span>mod<span class="p">)[</span><span class="m">-1</span><span class="p">],</span>    </div><div class='line' id='LC605'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;eval<span class="p">,</span> envir<span class="o">=</span><span class="kc">NULL</span><span class="p">,</span> enclos<span class="o">=</span><span class="kc">NULL</span><span class="p">))</span></div><div class='line' id='LC606'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">return</span><span class="p">(</span> list<span class="p">(</span>label<span class="o">=</span>label<span class="p">)</span> <span class="p">)</span></div><div class='line' id='LC607'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span> <span class="kr">else</span> <span class="kr">if</span><span class="p">(</span>mod<span class="p">[[</span><span class="m">1</span>L<span class="p">]]</span> <span class="o">==</span> <span class="s">&quot;label&quot;</span><span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC608'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;label <span class="o">&lt;-</span> unlist<span class="p">(</span>lapply<span class="p">(</span>as.list<span class="p">(</span>mod<span class="p">)[</span><span class="m">-1</span><span class="p">],</span></div><div class='line' id='LC609'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;eval<span class="p">,</span> envir<span class="o">=</span><span class="kc">NULL</span><span class="p">,</span> enclos<span class="o">=</span><span class="kc">NULL</span><span class="p">))</span></div><div class='line' id='LC610'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;label<span class="p">[</span>is.na<span class="p">(</span>label<span class="p">)]</span> <span class="o">&lt;-</span> <span class="s">&quot;&quot;</span> <span class="c1"># catch &#39;NA&#39; elements in a label</span></div><div class='line' id='LC611'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">return</span><span class="p">(</span> list<span class="p">(</span>label<span class="o">=</span>label<span class="p">)</span> <span class="p">)</span></div><div class='line' id='LC612'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span> <span class="kr">else</span> <span class="kr">if</span><span class="p">(</span>mod<span class="p">[[</span><span class="m">1</span>L<span class="p">]]</span> <span class="o">==</span> <span class="s">&quot;c&quot;</span><span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC613'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># vector: we allow numeric and character only!</span></div><div class='line' id='LC614'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;cof <span class="o">&lt;-</span> unlist<span class="p">(</span>lapply<span class="p">(</span>as.list<span class="p">(</span>mod<span class="p">)[</span><span class="m">-1</span><span class="p">],</span>    </div><div class='line' id='LC615'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;eval<span class="p">,</span> envir<span class="o">=</span><span class="kc">NULL</span><span class="p">,</span> enclos<span class="o">=</span><span class="kc">NULL</span><span class="p">))</span></div><div class='line' id='LC616'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>is.numeric<span class="p">(</span>cof<span class="p">))</span> </div><div class='line' id='LC617'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">return</span><span class="p">(</span> list<span class="p">(</span>fixed<span class="o">=</span>cof<span class="p">)</span> <span class="p">)</span></div><div class='line' id='LC618'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">else</span> <span class="kr">if</span><span class="p">(</span>is.character<span class="p">(</span>cof<span class="p">))</span> <span class="p">{</span></div><div class='line' id='LC619'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;cof<span class="p">[</span>is.na<span class="p">(</span>cof<span class="p">)]</span> <span class="o">&lt;-</span> <span class="s">&quot;&quot;</span> <span class="c1"># catch &#39;NA&#39; elements in a label</span></div><div class='line' id='LC620'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">return</span><span class="p">(</span> list<span class="p">(</span>label<span class="o">=</span>cof<span class="p">)</span> <span class="p">)</span></div><div class='line' id='LC621'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span> <span class="kr">else</span> <span class="p">{</span></div><div class='line' id='LC622'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;stop<span class="p">(</span><span class="s">&quot;lavaan ERROR: can not parse modifier:&quot;</span><span class="p">,</span> mod<span class="p">,</span> <span class="s">&quot;\n&quot;</span><span class="p">)</span></div><div class='line' id='LC623'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC624'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span> <span class="kr">else</span> <span class="p">{</span></div><div class='line' id='LC625'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># unknown expression</span></div><div class='line' id='LC626'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># as a final attempt, we will evaluate it and coerce it</span></div><div class='line' id='LC627'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># to either a numeric or character (vector)</span></div><div class='line' id='LC628'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;cof <span class="o">&lt;-</span> try<span class="p">(</span> eval<span class="p">(</span>mod<span class="p">,</span> envir<span class="o">=</span><span class="kc">NULL</span><span class="p">,</span> enclos<span class="o">=</span><span class="kc">NULL</span><span class="p">),</span> silent<span class="o">=</span><span class="kc">TRUE</span><span class="p">)</span></div><div class='line' id='LC629'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>is.numeric<span class="p">(</span>cof<span class="p">))</span></div><div class='line' id='LC630'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">return</span><span class="p">(</span> list<span class="p">(</span>fixed<span class="o">=</span>cof<span class="p">)</span> <span class="p">)</span></div><div class='line' id='LC631'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">else</span> <span class="kr">if</span><span class="p">(</span>is.character<span class="p">(</span>cof<span class="p">))</span></div><div class='line' id='LC632'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">return</span><span class="p">(</span> list<span class="p">(</span>label<span class="o">=</span>cof<span class="p">)</span> <span class="p">)</span></div><div class='line' id='LC633'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">else</span> <span class="p">{</span></div><div class='line' id='LC634'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;stop<span class="p">(</span><span class="s">&quot;lavaan ERROR: can not parse modifier:&quot;</span><span class="p">,</span> mod<span class="p">,</span> <span class="s">&quot;\n&quot;</span><span class="p">)</span></div><div class='line' id='LC635'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC636'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC637'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC638'><br/></div><div class='line' id='LC639'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># fill in rhs list</span></div><div class='line' id='LC640'>&nbsp;&nbsp;&nbsp;&nbsp;out <span class="o">&lt;-</span> list<span class="p">()</span></div><div class='line' id='LC641'>&nbsp;&nbsp;&nbsp;&nbsp;repeat <span class="p">{</span></div><div class='line' id='LC642'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>length<span class="p">(</span>rhs<span class="p">)</span> <span class="o">==</span> <span class="m">1</span>L<span class="p">)</span> <span class="p">{</span> <span class="c1"># last one and only a single element</span></div><div class='line' id='LC643'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;out <span class="o">&lt;-</span> c<span class="p">(</span>vector<span class="p">(</span><span class="s">&quot;list&quot;</span><span class="p">,</span> <span class="m">1</span>L<span class="p">),</span> out<span class="p">)</span></div><div class='line' id='LC644'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;NAME <span class="o">&lt;-</span> all.vars<span class="p">(</span>rhs<span class="p">)</span></div><div class='line' id='LC645'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>length<span class="p">(</span>NAME<span class="p">)</span> <span class="o">&gt;</span> <span class="m">0</span>L<span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC646'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;names<span class="p">(</span>out<span class="p">)[</span><span class="m">1</span>L<span class="p">]</span> <span class="o">&lt;-</span> NAME</div><div class='line' id='LC647'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span> <span class="kr">else</span> <span class="p">{</span> <span class="c1"># intercept</span></div><div class='line' id='LC648'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;names<span class="p">(</span>out<span class="p">)[</span><span class="m">1</span>L<span class="p">]</span> <span class="o">&lt;-</span> <span class="s">&quot;intercept&quot;</span></div><div class='line' id='LC649'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>as.character<span class="p">(</span>rhs<span class="p">)</span> <span class="o">==</span> <span class="s">&quot;0&quot;</span><span class="p">)</span></div><div class='line' id='LC650'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;out<span class="p">[[</span><span class="m">1</span>L<span class="p">]]$</span>fixed <span class="o">&lt;-</span> <span class="m">0</span></div><div class='line' id='LC651'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC652'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">break</span></div><div class='line' id='LC653'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span> <span class="kr">else</span> <span class="kr">if</span><span class="p">(</span>rhs<span class="p">[[</span><span class="m">1</span>L<span class="p">]]</span> <span class="o">==</span> <span class="s">&quot;*&quot;</span><span class="p">)</span> <span class="p">{</span> <span class="c1"># last one, but with modifier</span></div><div class='line' id='LC654'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;out <span class="o">&lt;-</span> c<span class="p">(</span>vector<span class="p">(</span><span class="s">&quot;list&quot;</span><span class="p">,</span> <span class="m">1</span>L<span class="p">),</span> out<span class="p">)</span></div><div class='line' id='LC655'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;NAME <span class="o">&lt;-</span> all.vars<span class="p">(</span>rhs<span class="p">[[</span><span class="m">3</span>L<span class="p">]])</span></div><div class='line' id='LC656'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>length<span class="p">(</span>NAME<span class="p">)</span> <span class="o">&gt;</span> <span class="m">0</span>L<span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC657'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;names<span class="p">(</span>out<span class="p">)[</span><span class="m">1</span>L<span class="p">]</span> <span class="o">&lt;-</span> NAME</div><div class='line' id='LC658'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span> <span class="kr">else</span> <span class="p">{</span> <span class="c1"># intercept</span></div><div class='line' id='LC659'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;names<span class="p">(</span>out<span class="p">)[</span><span class="m">1</span>L<span class="p">]</span> <span class="o">&lt;-</span> <span class="s">&quot;intercept&quot;</span></div><div class='line' id='LC660'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC661'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;i.var <span class="o">&lt;-</span> all.vars<span class="p">(</span>rhs<span class="p">[[</span><span class="m">2</span>L<span class="p">]],</span> unique<span class="o">=</span><span class="kc">FALSE</span><span class="p">)</span></div><div class='line' id='LC662'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>length<span class="p">(</span>i.var<span class="p">)</span> <span class="o">&gt;</span> <span class="m">0</span>L<span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC663'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># modifier are unquoted labels</span></div><div class='line' id='LC664'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;out<span class="p">[[</span><span class="m">1</span>L<span class="p">]]$</span>label <span class="o">&lt;-</span> i.var</div><div class='line' id='LC665'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span> <span class="kr">else</span> <span class="p">{</span></div><div class='line' id='LC666'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># modifer is something else</span></div><div class='line' id='LC667'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;out<span class="p">[[</span><span class="m">1</span>L<span class="p">]]</span> <span class="o">&lt;-</span> getModifier<span class="p">(</span>rhs<span class="p">[[</span><span class="m">2</span>L<span class="p">]])</span></div><div class='line' id='LC668'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC669'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">break</span></div><div class='line' id='LC670'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span> <span class="kr">else</span> <span class="kr">if</span><span class="p">(</span>rhs<span class="p">[[</span><span class="m">1</span>L<span class="p">]]</span> <span class="o">==</span> <span class="s">&quot;+&quot;</span><span class="p">)</span> <span class="p">{</span> <span class="c1"># not last one!</span></div><div class='line' id='LC671'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;i.var <span class="o">&lt;-</span> all.vars<span class="p">(</span>rhs<span class="p">[[</span><span class="m">3</span>L<span class="p">]],</span> unique<span class="o">=</span><span class="kc">FALSE</span><span class="p">)</span></div><div class='line' id='LC672'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;n.var <span class="o">&lt;-</span> length<span class="p">(</span>i.var<span class="p">)</span></div><div class='line' id='LC673'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;out <span class="o">&lt;-</span> c<span class="p">(</span>vector<span class="p">(</span><span class="s">&quot;list&quot;</span><span class="p">,</span> <span class="m">1</span>L<span class="p">),</span> out<span class="p">)</span></div><div class='line' id='LC674'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>length<span class="p">(</span>i.var<span class="p">)</span> <span class="o">&gt;</span> <span class="m">0</span>L<span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC675'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;names<span class="p">(</span>out<span class="p">)[</span><span class="m">1</span>L<span class="p">]</span> <span class="o">&lt;-</span> i.var<span class="p">[</span>n.var<span class="p">]</span></div><div class='line' id='LC676'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span> <span class="kr">else</span> <span class="p">{</span></div><div class='line' id='LC677'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;names<span class="p">(</span>out<span class="p">)[</span><span class="m">1</span>L<span class="p">]</span> <span class="o">&lt;-</span> <span class="s">&quot;intercept&quot;</span></div><div class='line' id='LC678'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC679'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>n.var <span class="o">&gt;</span> <span class="m">1</span>L<span class="p">)</span> <span class="p">{</span> </div><div class='line' id='LC680'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># modifier are unquoted labels</span></div><div class='line' id='LC681'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;out<span class="p">[[</span><span class="m">1</span>L<span class="p">]]$</span>label <span class="o">&lt;-</span> i.var<span class="p">[</span><span class="o">-</span>n.var<span class="p">]</span></div><div class='line' id='LC682'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span> <span class="kr">else</span> <span class="kr">if</span><span class="p">(</span>length<span class="p">(</span>rhs<span class="p">[[</span><span class="m">3</span><span class="p">]])</span> <span class="o">==</span> <span class="m">3</span>L<span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC683'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># modifiers!!</span></div><div class='line' id='LC684'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;out<span class="p">[[</span><span class="m">1</span>L<span class="p">]]</span> <span class="o">&lt;-</span> getModifier<span class="p">(</span>rhs<span class="p">[[</span><span class="m">3</span>L<span class="p">]][[</span><span class="m">2</span>L<span class="p">]])</span></div><div class='line' id='LC685'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC686'><br/></div><div class='line' id='LC687'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># next element</span></div><div class='line' id='LC688'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;rhs <span class="o">&lt;-</span> rhs<span class="p">[[</span><span class="m">2</span>L<span class="p">]]</span></div><div class='line' id='LC689'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span> <span class="kr">else</span> <span class="p">{</span></div><div class='line' id='LC690'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;stop<span class="p">(</span><span class="s">&quot;lavaan ERROR: I&#39;m confused parsing this line: &quot;</span><span class="p">,</span> rhs<span class="p">,</span> <span class="s">&quot;\n&quot;</span><span class="p">)</span></div><div class='line' id='LC691'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC692'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC693'><br/></div><div class='line' id='LC694'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1"># if multiple elements, check for duplicated elements and merge if found</span></div><div class='line' id='LC695'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">if</span><span class="p">(</span>length<span class="p">(</span>out<span class="p">)</span> <span class="o">&gt;</span> <span class="m">1</span>L<span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC696'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;rhs.names <span class="o">&lt;-</span> names<span class="p">(</span>out<span class="p">)</span></div><div class='line' id='LC697'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kr">while</span><span class="p">(</span> !is.na<span class="p">(</span>idx <span class="o">&lt;-</span> which<span class="p">(</span>duplicated<span class="p">(</span>rhs.names<span class="p">))[</span><span class="m">1</span>L<span class="p">])</span> <span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC698'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;dup.name <span class="o">&lt;-</span> rhs.names<span class="p">[</span> idx <span class="p">]</span></div><div class='line' id='LC699'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;orig.idx <span class="o">&lt;-</span> match<span class="p">(</span>dup.name<span class="p">,</span> rhs.names<span class="p">)</span></div><div class='line' id='LC700'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;out<span class="p">[[</span>orig.idx<span class="p">]]</span> <span class="o">&lt;-</span> c<span class="p">(</span> out<span class="p">[[</span>orig.idx<span class="p">]],</span> out<span class="p">[[</span>idx<span class="p">]]</span> <span class="p">)</span></div><div class='line' id='LC701'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;out <span class="o">&lt;-</span> out<span class="p">[</span><span class="o">-</span>idx<span class="p">]</span></div><div class='line' id='LC702'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;rhs.names <span class="o">&lt;-</span> names<span class="p">(</span>out<span class="p">)</span></div><div class='line' id='LC703'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC704'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC705'><br/></div><div class='line' id='LC706'>&nbsp;&nbsp;&nbsp;&nbsp;out</div><div class='line' id='LC707'><span class="p">}</span></div><div class='line' id='LC708'><br/></div></pre></div>
          </td>
        </tr>
      </table>
  </div>

          </div>
        </div>
      </div>
    </div>

  </div>

<div class="frame frame-loading large-loading-area" style="display:none;" data-tree-list-url="/yrosseel/lavaan/tree-list/79b51b6496090dca20de417c2b60c57d3f0c93d2" data-blob-url-prefix="/yrosseel/lavaan/blob/79b51b6496090dca20de417c2b60c57d3f0c93d2">
  <img src="https://a248.e.akamai.net/assets.github.com/images/spinners/octocat-spinner-64.gif?1334862346" height="64" width="64">
</div>

      </div>
      <div class="context-overlay"></div>
    </div>

      <div id="footer-push"></div><!-- hack for sticky footer -->
    </div><!-- end of wrapper - hack for sticky footer -->

      <!-- footer -->
      <div id="footer" >
        
  <div class="upper_footer">
     <div class="container clearfix">

       <!--[if IE]><h4 id="blacktocat_ie">GitHub Links</h4><![endif]-->
       <![if !IE]><h4 id="blacktocat">GitHub Links</h4><![endif]>

       <ul class="footer_nav">
         <h4>GitHub</h4>
         <li><a href="https://github.com/about">About</a></li>
         <li><a href="https://github.com/blog">Blog</a></li>
         <li><a href="https://github.com/features">Features</a></li>
         <li><a href="https://github.com/contact">Contact &amp; Support</a></li>
         <li><a href="https://github.com/training">Training</a></li>
         <li><a href="http://enterprise.github.com/">GitHub Enterprise</a></li>
         <li><a href="http://status.github.com/">Site Status</a></li>
       </ul>

       <ul class="footer_nav">
         <h4>Tools</h4>
         <li><a href="http://get.gaug.es/">Gauges: Analyze web traffic</a></li>
         <li><a href="http://speakerdeck.com">Speaker Deck: Presentations</a></li>
         <li><a href="https://gist.github.com">Gist: Code snippets</a></li>
         <li><a href="http://mac.github.com/">GitHub for Mac</a></li>
         <li><a href="http://windows.github.com/">GitHub for Windows</a></li>
         <li><a href="http://mobile.github.com/">Issues for iPhone</a></li>
         <li><a href="http://jobs.github.com/">Job Board</a></li>
       </ul>

       <ul class="footer_nav">
         <h4>Extras</h4>
         <li><a href="http://shop.github.com/">GitHub Shop</a></li>
         <li><a href="http://octodex.github.com/">The Octodex</a></li>
       </ul>

       <ul class="footer_nav">
         <h4>Documentation</h4>
         <li><a href="http://help.github.com/">GitHub Help</a></li>
         <li><a href="http://developer.github.com/">Developer API</a></li>
         <li><a href="http://github.github.com/github-flavored-markdown/">GitHub Flavored Markdown</a></li>
         <li><a href="http://pages.github.com/">GitHub Pages</a></li>
       </ul>

     </div><!-- /.site -->
  </div><!-- /.upper_footer -->

<div class="lower_footer">
  <div class="container clearfix">
    <!--[if IE]><div id="legal_ie"><![endif]-->
    <![if !IE]><div id="legal"><![endif]>
      <ul>
          <li><a href="https://github.com/site/terms">Terms of Service</a></li>
          <li><a href="https://github.com/site/privacy">Privacy</a></li>
          <li><a href="https://github.com/security">Security</a></li>
      </ul>

      <p>&copy; 2012 <span title="0.33890s from fe13.rs.github.com">GitHub</span> Inc. All rights reserved.</p>
    </div><!-- /#legal or /#legal_ie-->

      <div class="sponsor">
        <a href="http://www.rackspace.com" class="logo">
          <img alt="Dedicated Server" height="36" src="https://a248.e.akamai.net/assets.github.com/images/modules/footer/rackspaces_logo.png?1334862346" width="38" />
        </a>
        Powered by the <a href="http://www.rackspace.com ">Dedicated
        Servers</a> and<br/> <a href="http://www.rackspacecloud.com">Cloud
        Computing</a> of Rackspace Hosting<span>&reg;</span>
      </div>
  </div><!-- /.site -->
</div><!-- /.lower_footer -->

      </div><!-- /#footer -->

    

<div id="keyboard_shortcuts_pane" class="instapaper_ignore readability-extra" style="display:none">
  <h2>Keyboard Shortcuts <small><a href="#" class="js-see-all-keyboard-shortcuts">(see all)</a></small></h2>

  <div class="columns threecols">
    <div class="column first">
      <h3>Site wide shortcuts</h3>
      <dl class="keyboard-mappings">
        <dt>s</dt>
        <dd>Focus site search</dd>
      </dl>
      <dl class="keyboard-mappings">
        <dt>?</dt>
        <dd>Bring up this help dialog</dd>
      </dl>
    </div><!-- /.column.first -->

    <div class="column middle" style='display:none'>
      <h3>Commit list</h3>
      <dl class="keyboard-mappings">
        <dt>j</dt>
        <dd>Move selection down</dd>
      </dl>
      <dl class="keyboard-mappings">
        <dt>k</dt>
        <dd>Move selection up</dd>
      </dl>
      <dl class="keyboard-mappings">
        <dt>c <em>or</em> o <em>or</em> enter</dt>
        <dd>Open commit</dd>
      </dl>
      <dl class="keyboard-mappings">
        <dt>y</dt>
        <dd>Expand URL to its canonical form</dd>
      </dl>
    </div><!-- /.column.first -->

    <div class="column last" style='display:none'>
      <h3>Pull request list</h3>
      <dl class="keyboard-mappings">
        <dt>j</dt>
        <dd>Move selection down</dd>
      </dl>
      <dl class="keyboard-mappings">
        <dt>k</dt>
        <dd>Move selection up</dd>
      </dl>
      <dl class="keyboard-mappings">
        <dt>o <em>or</em> enter</dt>
        <dd>Open issue</dd>
      </dl>
      <dl class="keyboard-mappings">
        <dt><span class="platform-mac">⌘</span><span class="platform-other">ctrl</span> <em>+</em> enter</dt>
        <dd>Submit comment</dd>
      </dl>
    </div><!-- /.columns.last -->

  </div><!-- /.columns.equacols -->

  <div style='display:none'>
    <div class="rule"></div>

    <h3>Issues</h3>

    <div class="columns threecols">
      <div class="column first">
        <dl class="keyboard-mappings">
          <dt>j</dt>
          <dd>Move selection down</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>k</dt>
          <dd>Move selection up</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>x</dt>
          <dd>Toggle selection</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>o <em>or</em> enter</dt>
          <dd>Open issue</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt><span class="platform-mac">⌘</span><span class="platform-other">ctrl</span> <em>+</em> enter</dt>
          <dd>Submit comment</dd>
        </dl>
      </div><!-- /.column.first -->
      <div class="column last">
        <dl class="keyboard-mappings">
          <dt>c</dt>
          <dd>Create issue</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>l</dt>
          <dd>Create label</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>i</dt>
          <dd>Back to inbox</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>u</dt>
          <dd>Back to issues</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>/</dt>
          <dd>Focus issues search</dd>
        </dl>
      </div>
    </div>
  </div>

  <div style='display:none'>
    <div class="rule"></div>

    <h3>Issues Dashboard</h3>

    <div class="columns threecols">
      <div class="column first">
        <dl class="keyboard-mappings">
          <dt>j</dt>
          <dd>Move selection down</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>k</dt>
          <dd>Move selection up</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>o <em>or</em> enter</dt>
          <dd>Open issue</dd>
        </dl>
      </div><!-- /.column.first -->
    </div>
  </div>

  <div style='display:none'>
    <div class="rule"></div>

    <h3>Network Graph</h3>
    <div class="columns equacols">
      <div class="column first">
        <dl class="keyboard-mappings">
          <dt><span class="badmono">←</span> <em>or</em> h</dt>
          <dd>Scroll left</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt><span class="badmono">→</span> <em>or</em> l</dt>
          <dd>Scroll right</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt><span class="badmono">↑</span> <em>or</em> k</dt>
          <dd>Scroll up</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt><span class="badmono">↓</span> <em>or</em> j</dt>
          <dd>Scroll down</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>t</dt>
          <dd>Toggle visibility of head labels</dd>
        </dl>
      </div><!-- /.column.first -->
      <div class="column last">
        <dl class="keyboard-mappings">
          <dt>shift <span class="badmono">←</span> <em>or</em> shift h</dt>
          <dd>Scroll all the way left</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>shift <span class="badmono">→</span> <em>or</em> shift l</dt>
          <dd>Scroll all the way right</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>shift <span class="badmono">↑</span> <em>or</em> shift k</dt>
          <dd>Scroll all the way up</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>shift <span class="badmono">↓</span> <em>or</em> shift j</dt>
          <dd>Scroll all the way down</dd>
        </dl>
      </div><!-- /.column.last -->
    </div>
  </div>

  <div >
    <div class="rule"></div>
    <div class="columns threecols">
      <div class="column first" >
        <h3>Source Code Browsing</h3>
        <dl class="keyboard-mappings">
          <dt>t</dt>
          <dd>Activates the file finder</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>l</dt>
          <dd>Jump to line</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>w</dt>
          <dd>Switch branch/tag</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>y</dt>
          <dd>Expand URL to its canonical form</dd>
        </dl>
      </div>
    </div>
  </div>

  <div style='display:none'>
    <div class="rule"></div>
    <div class="columns threecols">
      <div class="column first">
        <h3>Browsing Commits</h3>
        <dl class="keyboard-mappings">
          <dt><span class="platform-mac">⌘</span><span class="platform-other">ctrl</span> <em>+</em> enter</dt>
          <dd>Submit comment</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>escape</dt>
          <dd>Close form</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>p</dt>
          <dd>Parent commit</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>o</dt>
          <dd>Other parent commit</dd>
        </dl>
      </div>
    </div>
  </div>
</div>

    <div id="markdown-help" class="instapaper_ignore readability-extra">
  <h2>Markdown Cheat Sheet</h2>

  <div class="cheatsheet-content">

  <div class="mod">
    <div class="col">
      <h3>Format Text</h3>
      <p>Headers</p>
      <pre>
# This is an &lt;h1&gt; tag
## This is an &lt;h2&gt; tag
###### This is an &lt;h6&gt; tag</pre>
     <p>Text styles</p>
     <pre>
*This text will be italic*
_This will also be italic_
**This text will be bold**
__This will also be bold__

*You **can** combine them*
</pre>
    </div>
    <div class="col">
      <h3>Lists</h3>
      <p>Unordered</p>
      <pre>
* Item 1
* Item 2
  * Item 2a
  * Item 2b</pre>
     <p>Ordered</p>
     <pre>
1. Item 1
2. Item 2
3. Item 3
   * Item 3a
   * Item 3b</pre>
    </div>
    <div class="col">
      <h3>Miscellaneous</h3>
      <p>Images</p>
      <pre>
![GitHub Logo](/images/logo.png)
Format: ![Alt Text](url)
</pre>
     <p>Links</p>
     <pre>
http://github.com - automatic!
[GitHub](http://github.com)</pre>
<p>Blockquotes</p>
     <pre>
As Kanye West said:

> We're living the future so
> the present is our past.
</pre>
    </div>
  </div>
  <div class="rule"></div>

  <h3>Code Examples in Markdown</h3>
  <div class="col">
      <p>Syntax highlighting with <a href="http://github.github.com/github-flavored-markdown/" title="GitHub Flavored Markdown" target="_blank">GFM</a></p>
      <pre>
```javascript
function fancyAlert(arg) {
  if(arg) {
    $.facebox({div:'#foo'})
  }
}
```</pre>
    </div>
    <div class="col">
      <p>Or, indent your code 4 spaces</p>
      <pre>
Here is a Python code example
without syntax highlighting:

    def foo:
      if not bar:
        return true</pre>
    </div>
    <div class="col">
      <p>Inline code for comments</p>
      <pre>
I think you should use an
`&lt;addr&gt;` element here instead.</pre>
    </div>
  </div>

  </div>
</div>


    <div class="ajax-error-message">
      <p><span class="mini-icon mini-icon-exclamation"></span> Something went wrong with that request. Please try again. <a href="javascript:;" class="ajax-error-dismiss">Dismiss</a></p>
    </div>

    <div id="logo-popup">
      <h2>Looking for the GitHub logo?</h2>
      <ul>
        <li>
          <h4>GitHub Logo</h4>
          <a href="http://github-media-downloads.s3.amazonaws.com/GitHub_Logos.zip"><img alt="Github_logo" src="https://a248.e.akamai.net/assets.github.com/images/modules/about_page/github_logo.png?1334862345" /></a>
          <a href="http://github-media-downloads.s3.amazonaws.com/GitHub_Logos.zip" class="minibutton btn-download download"><span class="icon"></span>Download</a>
        </li>
        <li>
          <h4>The Octocat</h4>
          <a href="http://github-media-downloads.s3.amazonaws.com/Octocats.zip"><img alt="Octocat" src="https://a248.e.akamai.net/assets.github.com/images/modules/about_page/octocat.png?1334862345" /></a>
          <a href="http://github-media-downloads.s3.amazonaws.com/Octocats.zip" class="minibutton btn-download download"><span class="icon"></span>Download</a>
        </li>
      </ul>
    </div>

    
    
    
    <span id='server_response_time' data-time='0.34185' data-host='fe13'></span>
  </body>
</html>

